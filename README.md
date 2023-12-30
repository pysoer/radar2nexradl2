# radar2nexradl2
将中国气象雷达标准格式转成美国nexradL2格式（.ar2v），可供GR2Analyst读取。

# 来源程序
来自一个南京大学的人写的，源代码运行有些许错误，本人只是做了修改。
[链接地址](https://git.nju.edu.cn/bofan/radar_read_write_plot)


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
