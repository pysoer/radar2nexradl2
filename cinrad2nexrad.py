# Edit: pysoer,QQ-group:480305660
# Script: Radar data format transfer: CINRAD format -> NEXRAD level II ar2v format
# Author: bofan@smail.nju.edu.cn (Thanks for any bug report)
# Version: 20230912 v1.0
# Path: D:\Seafile\Share\notebooks\20230601ReadFile\CINRAD2NEXRAD.py
import struct
import pyart.io.nexrad_level2 as ar2v
from cinradReader import cinradReader, data_type
from pyart.io.nexrad_level2 import _structure_size
import warnings


def cinrad2nexrad(inpath, outpath, site_name_4s=None):
    # 1. Read CINRAD radar product
    radar = cinradReader(inpath)
    # 2. Write NEXRAD level II product
    ## 2.1 Write volume header
    volume_header = []  # ar2v.volume_header
    volume_header = convert_to_dict(ar2v.VOLUME_HEADER)
    volume_header["tape"] = "AR2V0006."
    volume_header["extension"] = "001"
    tt = radar.info["task"]["scan_start_time"]
    volume_header["date"] = int(tt / 24 / 3600)
    volume_header["time"] = int((tt / 24 / 3600 - int(tt / 24 / 3600)) * 24 * 3600 * 1000)
    if site_name_4s is None:
        site_name_4s = str(radar.info["site"]["site_name"][0:4], encoding="ASCII")
    volume_header[
        "icao"
    ] = site_name_4s  # str(radar.info['site']['site_name'][0:4],encoding='ASCII')

    rec = []  # various dict of 'header','msg_header','VOL','ELV','RAD' and other fields
    list_index = [2, 3, 7, 10, 9]  # REF, VEL, ZDR, PHI, RHO

    print("rec-data-len:",len(radar.rec["data"]))
    for i in range(len(radar.rec["data"])):
        rec.append({})

        ## 2.2.2 Write data
        for j in ["VOL", "ELV", "RAD"]:
            field = j
            if field == "VOL":
                block = convert_to_dict(ar2v.VOLUME_DATA_BLOCK)
                block["block_type"] = "R"
                block["data_name"] = field
                block["lrtup"] = 44
                block["version_major"] = 1
                block["version_minor"] = 0
                block["lat"] = radar.info["site"]["Latitude"][0]
                block["lon"] = radar.info["site"]["Longitude"][0]
                block["height"] = radar.info["site"]["antenna_height"][0]
                block["feedhorn_height"] = radar.info["site"]["antenna_height"][0]  # ?
                block["refl_calib"] = radar.info["task"]["hori_cali"][
                    0
                ]  # ? e.g. KCLX: -44.6
                block["power_h"] = 269.2  # ? e.g. KCLX
                block["power_v"] = 266.7  # ? e.g. KCLX
                block["diff_refl_calib"] = radar.info["task"]["ZDR_cali"][
                    0
                ]  # e.g. KCLX: 0.35
                block["init_phase"] = (
                    radar.info["task"]["PHIDP_cali"] + 180
                )  # e.g. KCLX: 60
                block["vcp"] = 212
                block["spare"] = ""  # ?
            elif field == "ELV":
                block = convert_to_dict(ar2v.ELEVATION_DATA_BLOCK)
                block["block_type"] = "R"
                block["data_name"] = field
                block["lrtup"] = 12
                block["atmos"] = int(
                    radar.info["cut"][
                        radar.rec["header"][i]["elevation_number"][0] - 1
                    ]["atmos_loss"]
                )
                block["refl_calib"] = radar.info["task"]["hori_cali"][
                    0
                ]  # ? e.g. KCLX: -43.8
            elif field == "RAD":
                block = convert_to_dict(ar2v.RADIAL_DATA_BLOCK)
                block["block_type"] = "R"
                block["data_name"] = field
                block["lrtup"] = 28
                block["unambig_range"] = int(
                    1.9e8
                    / radar.info["cut"][
                        radar.rec["header"][i]["elevation_number"][0] - 1
                    ]["PRF1"]
                    / 1000
                )
                block["noise_h"] = radar.info["task"]["hori_noise"][0]
                block["noise_v"] = radar.info["task"]["vert_noise"][0]
                block["nyquist_vel"] = int(
                    radar.info["cut"][
                        radar.rec["header"][i]["elevation_number"][0] - 1
                    ]["nyquist_spd"]
                )
                block["spare"] = ""  # ?
            rec[i][field] = block

        for j in range(len(radar.rec["data"][i]["moment_header"])):
            if radar.rec["data"][i]["moment_header"][j]["data_type"] not in list_index:
                continue
            field = data_type[radar.rec["data"][i]["moment_header"][j]["data_type"][0]]
            block = convert_to_dict(ar2v.GENERIC_DATA_BLOCK)
            block["block_type"] = "D"
            block["data_name"] = field
            block["reserved"] = 0  # ?
            block["ngates"] = int(
                radar.rec["data"][i]["moment_header"][j]["block_length"]
                / radar.rec["data"][i]["moment_header"][j]["bin_length"]
            )
            block["first_gate"] = radar.info["cut"][
                radar.rec["header"][i]["elevation_number"][0] - 1
            ]["start_range"][0]
            if field == "VEL":
                block["gate_spacing"] = radar.info["cut"][
                    radar.rec["header"][i]["elevation_number"][0] - 1
                ]["dop_reso"][0]
            else:
                block["gate_spacing"] = radar.info["cut"][
                    radar.rec["header"][i]["elevation_number"][0] - 1
                ]["log_reso"][0]
            block["thresh"] = 100  # ?
            block["snr_thres"] = 16  # ?
            block["flags"] = 0  # ?
            block["word_size"] = (
                8 * radar.rec["data"][i]["moment_header"][j]["bin_length"][0]
            )
            block["scale"] = radar.rec["data"][i]["moment_header"][j]["scale"][0]
            block["offset"] = radar.rec["data"][i]["moment_header"][j]["offset"][0]

            rec[i][field] = block
            rec[i][field]["data"] = radar.rec["data"][i]["moment_data"][j]

        ## 2.2.1 Write data header
        msg_header = convert_to_dict(ar2v.MSG_31)
        msg_header["id"] = volume_header["icao"]
        msg_header["collect_ms"] = int(
            (
                radar.rec["header"][i]["seconds"] / 24 / 3600
                - int(radar.rec["header"][i]["seconds"] / 24 / 3600)
            )
            * 24
            * 3600
            * 1000
            + radar.rec["header"][i]["microseconds"]
        )
        msg_header["collect_date"] = int(radar.rec["header"][i]["seconds"] / 24 / 3600)
        msg_header["azimuth_number"] = radar.rec["header"][i]["radial_number"][0]
        msg_header["azimuth_angle"] = radar.rec["header"][i]["azimuth"][0]
        msg_header["compress_flag"] = 0
        msg_header["spare_0"] = 0  # ?
        msg_header["radial_length"] = []
        msg_header["azimuth_resolution"] = 0  # ?
        msg_header["radial_spacing"] = radar.rec["header"][i]["radial_state"][
            0
        ]  # 0-4 the same, >4 wrong
        msg_header["elevation_number"] = radar.rec["header"][i]["elevation_number"][0]
        msg_header["cut_sector"] = 1
        msg_header["elevation_angle"] = radar.rec["header"][i]["elevation"][0]
        msg_header["radial_blanking"] = 0
        msg_header["azimuth_mode"] = 0
        msg_header["block_count"] = 0
        msg_header["block_pointer_1"] = 0
        msg_header["block_pointer_2"] = 0
        msg_header["block_pointer_3"] = 0
        msg_header["block_pointer_4"] = 0
        msg_header["block_pointer_5"] = 0
        msg_header["block_pointer_6"] = 0
        msg_header["block_pointer_7"] = 0
        msg_header["block_pointer_8"] = 0
        msg_header["block_pointer_9"] = 0
        msg_header["block_pointer_10"] = 0
        k = 1
        fields = sorted(
            list(rec[i].keys()),
            key=(["VOL", "ELV", "RAD"] + [data_type[ii] for ii in list_index]).index,
        )
        for field in fields:
            if k > 9:
                print("Too many fields(>10) to record!")
                break
            if field == "VOL" and k == 1:
                msg_header["block_pointer_1"] = _structure_size(ar2v.MSG_31)
            elif field == "ELV" and k == 2:
                msg_header["block_pointer_2"] = msg_header[
                    "block_pointer_1"
                ] + _structure_size(ar2v.VOLUME_DATA_BLOCK)
            elif field == "RAD" and k == 3:
                msg_header["block_pointer_3"] = msg_header[
                    "block_pointer_2"
                ] + _structure_size(ar2v.ELEVATION_DATA_BLOCK)
            elif k == 4:
                msg_header["block_pointer_4"] = msg_header[
                    "block_pointer_3"
                ] + _structure_size(ar2v.RADIAL_DATA_BLOCK)
            else:
                msg_header["block_pointer_" + str(k)] = (
                    msg_header["block_pointer_" + str(k - 1)]
                    + _structure_size(ar2v.GENERIC_DATA_BLOCK)
                    + int(rec[i][field]["word_size"] / 8 * rec[i][field]["ngates"])
                )
            k = k + 1
        msg_header["block_count"] = k - 1
        msg_header["radial_length"] = (
            msg_header["block_pointer_" + str(k - 1)]
            + _structure_size(ar2v.GENERIC_DATA_BLOCK)
            + int(rec[i][fields[-1]]["word_size"] / 8 * rec[i][fields[-1]]["ngates"])
        )
        rec[i]["msg_header"] = msg_header

        ## 2.2 Write message header
        header = convert_to_dict(ar2v.MSG_HEADER)
        header["size"] = int(
            (msg_header["radial_length"] + _structure_size(ar2v.MSG_HEADER)) / 2
        )  # ?
        header["channels"] = 0  # ?
        header["type"] = 31
        header["seq_id"] = 0  # ?
        header["date"] = 0  # ?
        header["ms"] = 0  # ?
        header["segments"] = 0  # ?
        header["seg_num"] = 0  # ?
        rec[i]["header"] = header

    # print(
    #     volume_header["icao"],
    #     volume_header["icao"],
    #     radar.info["site"]["Latitude"],
    #     radar.info["site"]["Longitude"],
    #     radar.info["site"]["antenna_height"],
    #     13,
    #     "CN",
    #     site_name_4s,
    #     "\n",
    # )

    # 3.Write into NEXRAD ar2v binary file

    with open(outpath, "wb") as f:
        # write Volume Header Record
        structure = ar2v.VOLUME_HEADER
        args = dic2vars(volume_header)
        lst = _pack_structure(structure, args)
        f.write(lst)

        # write Local Data Manager (LDM) Compressed Record
        ## Message Data
        for i in range(len(rec)):
            f.write(b"\x00" * 12)
            ### MSG Header
            structure = ar2v.MSG_HEADER
            args = dic2vars(rec[i]["header"])
            lst = _pack_structure(structure, args)
            f.write(lst)
            ### MSG 31 Header
            structure = ar2v.MSG_31
            args = dic2vars(rec[i]["msg_header"])
            lst = _pack_structure(structure, args)
            # f.write(lst3[:-4]) # Block_Pointer_10 wrong to RVOL

            list1 = [
                "VOL",
                "ELV",
                "RAD",
                "REF",
                "VEL",
                "SW",
                "ZDR",
                "PHI",
                "RHO",
                "CFP",
            ]
            set2 = rec[i].keys()
            list3 = list(list1 & set2)
            list3 = sorted(list3, key=list1.index)

            for j in list3:  # the sequence of fields in this MSG 31
                #### MSG 31 Data Block
                if j == "VOL":
                    structure = ar2v.VOLUME_DATA_BLOCK
                    args = dic2vars(rec[i][j])
                    new_lst = _pack_structure(structure, args)
                    lst = lst_plus(
                        lst,
                        new_lst,
                        rec[i]["msg_header"][
                            "block_pointer_" + str(list3.index(j) + 1)
                        ],
                    )
                elif j == "ELV":
                    structure = ar2v.ELEVATION_DATA_BLOCK
                    args = dic2vars(rec[i][j])
                    new_lst = _pack_structure(structure, args)
                    lst = lst_plus(
                        lst,
                        new_lst,
                        rec[i]["msg_header"][
                            "block_pointer_" + str(list3.index(j) + 1)
                        ],
                    )
                elif j == "RAD":
                    structure = ar2v.RADIAL_DATA_BLOCK
                    args = dic2vars(rec[i][j])
                    new_lst = _pack_structure(structure, args)
                    lst = lst_plus(
                        lst,
                        new_lst,
                        rec[i]["msg_header"][
                            "block_pointer_" + str(list3.index(j) + 1)
                        ],
                    )
                else:
                    structure = ar2v.GENERIC_DATA_BLOCK
                    tmp_rec = rec[i][j].copy()
                    tmp_rec.pop("data")
                    args = dic2vars(tmp_rec)
                    new_lst = _pack_structure(structure, args)
                    lst = lst_plus(
                        lst,
                        new_lst,
                        rec[i]["msg_header"][
                            "block_pointer_" + str(list3.index(j) + 1)
                        ],
                    )
                    ##### MSG 31 Generic Data
                    tmp_data = rec[i][j]["data"]
                    l, ngates = len(tmp_data), tmp_rec["ngates"]
                    if tmp_rec["word_size"] == 16:
                        tmp_format = ">" + "H" * ngates
                        lst = lst + struct.pack(tmp_format, *tmp_data)
                    elif tmp_rec["word_size"] == 8:
                        tmp_format = ">" + "B" * ngates
                        lst = lst + struct.pack(tmp_format, *tmp_data)
                    else:
                        warnings.warn(
                            'Unsupported bit size: %s. Returning array dtype "B"',
                            tmp_rec["word_size"],
                        )
            f.write(lst)
        f.closed
    return


def convert_to_dict(tuple_list):
    # Create a dictionary using the dict() constructor and a list comprehension
    dictionary = dict((key, value) for key, value in tuple_list)
    # Return the completed dictionary
    return dictionary


def dic2vars(dic):
    # print([dic[i] for i in dic.keys()])
    args = []
    for i in dic.keys():
        if isinstance(dic[i], str):
            args.append(bytes(dic[i], encoding="ASCII"))
        else:
            args.append(dic[i])
    # print(args)
    return args


def _pack_structure(structure, args):
    fmt = ">" + "".join(i[1] for i in structure)  # NEXRAD is big-endian
    lst = struct.pack(fmt, *args)
    return lst


def lst_plus(lst, new_lst, pointer):
    if len(lst) < pointer:
        lst = lst + b"\00" * (pointer - len(lst))
        lst = lst + new_lst
    elif len(lst) >= pointer:
        lst = lst[:pointer]
        lst = lst + new_lst
    return lst


if __name__ == "__main__":
    indir = (
        "D:/Seafile/Share/notebooks/20230523DisplayGUI/dual_pol_radar_plot/testdata/"
    )
    fname = "NJU.20220729.105032.AR2.bz2"
    inpath = indir + fname
    outdir = "./"
    outpath = outdir + fname + ".ar2v"
    cinrad2nexrad(inpath, outpath, "NJUC")
