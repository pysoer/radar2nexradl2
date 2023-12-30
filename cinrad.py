# Edit: pysoer,QQ-group:480305660
# Module: Read CINRAD radar binary file(.bin/.bin.bz2)
# Author: bofan@smail.nju.edu.cn(Thanks for any bug report)
# Version: 20230910 alpha, D:\Seafile\Share\notebooks\20230601ReadFile\CINRAD_test_read.ipynb
#          20230911 beta, D:\Seafile\Share\notebooks\20230601ReadFile\CINRAD_read.py

import numpy as np
import struct
import warnings
import re


class cinrad(object):
    """
    Class for accessing all info and data in a CINRAD (WSR-98D) binary file.

    CINRAD files [1]_, also know as WSR-98D weather radar base data.
    This class supports encoding the whole format of a CINRAD file.
    Files with uncompressed messages and compressed messages are supported.

    Parameters
    ----------
    inpath : str
        Path of CINRAD file to read.

    Attributes
    ----------
    info: list
        Comman blocks(generic, site, task and cuts) in the file.
    rec : list
        Radial blocks, including header and data(moment_header, moment_data), in the file.

    References
    ----------
    .. [1] https://www.cma.gov.cn/zfxxgk/gknr/flfgbz/bz/202307/t20230712_5642881.html
    .. [2] pyart.io.nexrad_level2
    .. [3] cinrad.io.StandardData

    """

    dtype_corr = {
        1: "TREF",
        2: "REF",
        3: "VEL",
        4: "SW",
        5: "SQI",
        6: "CPA",
        7: "ZDR",
        8: "LDR",
        9: "RHO",
        10: "PHI",
        11: "KDP",
        12: "CP",
        14: "HCL",
        15: "CF",
        16: "SNRH",
        17: "SNRV",
        19: "POTS",
        21: "COP",
        26: "VELSZ",
        27: "DR",
        32: "Zc",
        33: "Vc",
        34: "Wc",
        35: "ZDRc",
    }

    def __init__(self, inpath):
        with _prepare_file(inpath) as self.f:
            self._parse()

    def _parse(self):
        header = self._unpack_from_buf(self.f, GENERIC_HEADER)
        site = self._unpack_from_buf(self.f, SITE_CONFIG)
        task = self._unpack_from_buf(self.f, TASK_CONFIG)
        cut = []
        for i in range(task["cut_number"][0]):
            cut.append(self._unpack_from_buf(self.f, CUT_CONFIG))
        info = {"header": header, "site": site, "task": task, "cut": cut}
        rec = {"header": [], "data": []}
        i=0
        while 1:
            try:
                header_bytes = self.f.read(64)
                if not header_bytes:
                    # Fix for single-tilt file
                    break
                radial_header = np.frombuffer(header_bytes, RADIAL_HEADER)
                rec["header"].append(radial_header)
                rec["data"].append({"moment_header": [], "moment_data": []})
                for _ in range(radial_header["moment_number"][0]):
                    moment_header = self._unpack_from_buf(self.f, MOMENT_HEADER)
                    if moment_header["block_length"] == 0:
                        # Some ill-formed files
                        continue
                    blo_len = int(moment_header["block_length"][0])
                    dtype_code = moment_header["data_type"][0]
                    dtype = self.dtype_corr.get(dtype_code, None)
                    bin_len = moment_header["bin_length"][0]
                    moment_data = np.frombuffer(
                        self.f.read(blo_len), "u{}".format(bin_len)
                    )
                    if not dtype:
                        warnings.warn(
                            'Unsupported bit size: {}. Returning array dtype "B"'.format(
                                bin_len
                            )
                        )
                    rec["data"][i]["moment_header"].append(moment_header)
                    rec["data"][i]["moment_data"].append(moment_data)
                i+=1
            except:
                break
        # print('Reading: %s'%inpath)
        print("First radial num: %d \nLast radial num: %d" % (0, i - 1))
        print('All info & data in dict: "info" and "rec"\n')

        # raw
        self.info = info
        self.rec = rec

        # spherical
        self.data, self.gate, self.azimuth, self.elevation = self.get_data()

        # cartesian
        self.spherical2cartesian = spherical2cartesian

    def get_data(self):
        """
        input: cinrad object
        output: data, gate, azimuth, elevation
        """
        import numpy as np

        nrays = len(self.rec["data"])
        nswps = self.info["task"]["cut_number"][0]
        self.nswps = nswps
        print("Number of rays: %d" % nrays)
        print("Number of sweeps: %d" % nswps)
        # data={'elevation':[],'azimuth':[]}
        rays_each_swp = [[] for _ in range(nswps)]
        for nray in range(nrays):
            swp = self.rec["header"][nray]["elevation_number"][0]
            rays_each_swp[swp - 1].append(nray)

        data = [[] for _ in range(nswps)]
        gate = [[] for _ in range(nswps)]
        elevation = [[] for _ in range(nswps)]
        azimuth = [[] for _ in range(nswps)]

        for nswp in range(nswps):
            field_list = []
            for nray in rays_each_swp[nswp]:
                ele = self.rec["header"][nray]["elevation"][0]
                azi = self.rec["header"][nray]["azimuth"][0]
                # print('nswp=%d nray=%d'%(nswp,nray),'ele=%.1f azi=%.1f'%(ele,azi))
                elevation[nswp].append(ele)
                azimuth[nswp].append(azi)
                if nray == rays_each_swp[nswp][0]:
                    data[nswp] = {}
                    gate[nswp] = {}
                ll_header = self.rec["data"][nray]["moment_header"]
                for nvar in range(len(ll_header)):
                    field = data_type[
                        self.rec["data"][nray]["moment_header"][nvar]["data_type"][0]
                    ]
                    tmp_data = self.rec["data"][nray]["moment_data"][nvar].copy()
                    scale = self.rec["data"][nray]["moment_header"][nvar]["scale"][0]
                    offset = self.rec["data"][nray]["moment_header"][nvar]["offset"][0]
                    mask = tmp_data < 5  # int 0-4 with special meanings
                    tmp_data = np.ma.array((tmp_data - offset) / scale, mask=mask)
                    if (
                        self.info["cut"][nswp]["log_reso"]
                        != self.info["cut"][nswp]["log_reso"]
                    ):
                        print("Error in range reso!")
                    ngates = len(self.rec["data"][nray]["moment_data"][nvar])
                    tmp_gate = (
                        np.array(range(ngates)) * self.info["cut"][nswp]["log_reso"]
                        + self.info["cut"][nswp]["start_range"]
                    )
                    if field not in field_list:
                        gate[nswp][field] = []
                        data[nswp][field] = []
                        field_list.append(field)
                    data[nswp][field].append(tmp_data)
                    if len(tmp_gate) > len(gate[nswp][field]):
                        gate[nswp][field] = tmp_gate
                # print(field_list)
            for field in field_list:
                data[nswp][field] = np.ma.vstack(data[nswp][field])
        return data, gate, azimuth, elevation

    def get_cut(self, field="REF", nswp=0):  # sp for class cinrad
        """
        input: field name, cut number
        output: range, azimuth, elevation, field data
        """
        ran = self.gate[nswp][field]
        azi = self.azimuth[nswp]
        ele = self.elevation[nswp]
        field_data = self.data[nswp][field]
        return ran, azi, ele, field_data
    
    def _unpack_from_buf(self,a,structure):
        """Unpack a structure from a buffer."""
        # size = _structure_size(structure)
        # return _unpack_structure(buf[pos:pos + size], structure)
        ndtype = np.dtype(structure)
        bsize = self._structure_size(structure)
        return np.frombuffer(self.f.read(bsize), dtype=ndtype)
    
    def _structure_size(self,structure):
        """Find the size of a structure in bytes."""
        # return struct.calcsize(''.join([i[1] for i in structure]))
        total = 0
        for _, dtype in structure:
            match = re.search(r"\d+", dtype)
            if match:
                total += int(match.group())
        return total

def get_cut(data, gate, elevation, azimuth, field="REF", nswp=0):
    """
    input: data, gate, elevation, azimuth, field name, cut number
    output: range, azimuth, elevation, field data
    """
    ran = gate[nswp][field]
    azi = azimuth[nswp]
    ele = elevation[nswp]
    field_data = data[nswp][field]
    return ran, azi, ele, field_data


def spherical2cartesian(ran, azi, ele, field_data=None):
    """
    output: range, azimuth, elevation, field data
    output: xx(km), yy(km), cc
    """
    from pyart.core.transforms import antenna_vectors_to_cartesian

    x, y, z = antenna_vectors_to_cartesian(ran, azi, ele)
    xx, yy, cc = x / 1000, y / 1000, field_data
    return xx, yy, cc


def _prepare_file(inpath):
    # function from cinrad.io.base

    import bz2
    import gzip

    if hasattr(inpath, "read"):
        return inpath
    f = open(inpath, "rb")
    magic = f.read(3)
    f.close()
    if magic.startswith(b"\x1f\x8b"):
        return gzip.GzipFile(inpath, "rb")
    if magic.startswith(b"BZh"):
        return bz2.BZ2File(inpath, "rb")
    return open(inpath, "rb")


data_type = {
    1: "dBT",
    2: "REF",  #'dBZ',
    3: "VEL",  #'V',
    4: "W",
    5: "SQI",
    6: "CPA",
    7: "ZDR",
    8: "LDR",
    9: "RHO",
    10: "PHI",
    11: "KDP",
    12: "CP",
    14: "HCL",
    15: "CF",
    16: "SNRH",
    17: "SNRV",
    19: "POTS",
    21: "COP",
    26: "VELSZ",
    27: "DR",
    32: "Zc",
    33: "Vc",
    34: "Wc",
    35: "ZDRc",
}
GENERIC_HEADER = [
    ("magic_number", "i4"),
    ("major_version", "i2"),
    ("minor_version", "i2"),
    ("generic_type", "i4"),
    ("product_type", "i4"),
    ("res1", "16c"),
]

SITE_CONFIG = [
    ("site_code", "S8"),
    ("site_name", "S32"),
    ("Latitude", "f4"),
    ("Longitude", "f4"),
    ("antenna_height", "i4"),
    ("ground_height", "i4"),
    ("frequency", "f4"),
    ("beam_width_hori", "f4"),
    ("beam_width_vert", "f4"),
    ("RDA_version", "i4"),
    ("radar_type", "i2"),
    ("antenna_gain", "i2"),
    ("trans_loss", "i2"),
    ("recv_loss", "i2"),
    ("other_loss", "i2"),
    ("res2", "46c"),
]

TASK_CONFIG = [
    ("task_name", "S32"),
    ("task_dsc", "S128"),
    ("polar_type", "i4"),
    ("scan_type", "i4"),
    ("pulse_width", "i4"),
    ("scan_start_time", "i4"),
    ("cut_number", "i4"),
    ("hori_noise", "f4"),
    ("vert_noise", "f4"),
    ("hori_cali", "f4"),
    ("vert_cali", "f4"),
    ("hori_tmp", "f4"),
    ("vert_tmp", "f4"),
    ("ZDR_cali", "f4"),
    ("PHIDP_cali", "f4"),
    ("LDR_cali", "f4"),
    ("res3", "40c"),
]

CUT_CONFIG = [
    ("process_mode", "i4"),
    ("wave_form", "i4"),
    ("PRF1", "f4"),
    ("PRF2", "f4"),
    ("dealias_mode", "i4"),
    ("azimuth", "f4"),
    ("elev", "f4"),
    ("start_angle", "f4"),
    ("end_angle", "f4"),
    ("angular_reso", "f4"),
    ("scan_spd", "f4"),
    ("log_reso", "i4"),
    ("dop_reso", "i4"),
    ("max_range1", "i4"),
    ("max_range2", "i4"),
    ("start_range", "i4"),
    ("sample1", "i4"),
    ("sample2", "i4"),
    ("phase_mode", "i4"),
    ("atmos_loss", "f4"),
    ("nyquist_spd", "f4"),
    ("moments_mask", "i8"),
    ("moments_size_mask", "i8"),
    ("misc_filter_mask", "i4"),
    ("SQI_thres", "f4"),
    ("SIG_thres", "f4"),
    ("CSR_thres", "f4"),
    ("LOG_thres", "f4"),
    ("CPA_thres", "f4"),
    ("PMI_thres", "f4"),
    ("DPLOG_thres", "f4"),
    ("res_thres", "4V"),
    ("dBT_mask", "i4"),
    ("dBZ_mask", "i4"),
    ("vel_mask", "i4"),
    ("sw_mask", "i4"),
    ("DP_mask", "i4"),
    ("res_mask", "12V"),
    ("scan_sync", "i4"),
    ("direction", "i4"),
    ("ground_clutter_classifier_type", "i2"),
    ("ground_clutter_filter_type", "i2"),
    ("ground_clutter_filter_notch_width", "i2"),
    ("ground_clutter_filter_window", "i2"),
    ("res4", "72V"),
]

RADIAL_HEADER = [
    ("radial_state", "i4"),
    ("spot_blank", "i4"),
    ("seq_number", "i4"),
    ("radial_number", "i4"),
    ("elevation_number", "i4"),
    ("azimuth", "f4"),
    ("elevation", "f4"),
    ("seconds", "i4"),
    ("microseconds", "i4"),
    ("data_length", "i4"),
    ("moment_number", "i4"),
    ("res5", "i2"),
    ("hori_est_noise", "i2"),
    ("vert_est_noise", "i2"),
    ("zip_type", "c"),
    ("res6", "13c"),
]

MOMENT_HEADER = [
    ("data_type", "i4"),
    ("scale", "i4"),
    ("offset", "i4"),
    ("bin_length", "i2"),
    ("flags", "i2"),
    ("block_length", "i4"),
    ("res", "12c"),
]
