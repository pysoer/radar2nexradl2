# radar2nexradl2
将中国气象雷达标准格式转成美国nexradL2格式（.ar2v），可供GR2Analyst读取。

# 来源程序
来自一个南京大学的人写的，源代码运行有些许错误，本人只是做了修改。
[链接地址](https://git.nju.edu.cn/bofan/radar_read_write_plot)

## 使用方法
```python
from cinrad2nexrad import cinrad2nexrad
import os
nFiles = "Z_RADR_I_Z9999_20230504102859_O_DOR_SAD_CAP_FMT.bin.bz2"
basename = nFiles+".ar2v"
cinrad2nexrad(nFiles,basename,"YYLD") 
# YYLD为雷达站名，只能是四个大写字母，要和grlevel2.cfg里面的一致
```

## 转完了后可以用pyart验证一下
```python
import pyart
nFiles = "Z_RADR_I_Z9999_20230504102859_O_DOR_SAD_CAP_FMT.bin.bz2.ar2v"
radar = pyart.io.read_nexrad_archive(nFiles)
# print(radar.scan_info())
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure(figsize=(6, 5))

# plot super resolution reflectivity
ax = fig.add_subplot(111)
display.plot(
    "reflectivity",
    0,
    title="NEXRAD Reflectivity",
    vmin=-32,
    vmax=64,
    colorbar_label="",
    ax=ax,
)
display.plot_range_ring(radar.range["data"][-1] / 1000.0, ax=ax)
display.set_limits(xlim=(-500, 500), ylim=(-500, 500), ax=ax)
plt.show()
```


# 其他说明
Edit: pysoer,QQ-group:480305660  

Script: radar read/write/plot python toolbox
Author: bofan@smail.nju.edu.cn
Version: 20231026 v1.0

-What can it do?
1. read CINRAD/NEXRAD II raw data
2. convert CINRAD format to NEXRAD format
3. plot radar map

-How to use it?
all modules is <class> or <function> of detailed annotation

-PS
Thanks for any bug report! 

-Reference  

[1] https://www.cma.gov.cn/zfxxgk/gknr/flfgbz/bz/202307/t20230712_5642881.html  
[2] pyart.io.nexrad_level2  
[3] cinrad.io.StandardData  
[4] https://www.roc.noaa.gov/wsr88d/PublicDocs/ICDs/2620010G.pdf  
[5] https://www.roc.noaa.gov/wsr88d/PublicDocs/ICDs/2620002U.pdf  
[6] thanks for PKU CLY, CQX. read_ar2v_radar.m  

