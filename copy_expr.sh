#!/bin/sh

# dataTag
dataTag=$1

echo download data ${dataTag} ...

# set directory
# branch_dir="/g3/wanghao/kezuo/xhj/GRAPES_GFS3.2/GRAPES_GFS3.2_fix_autobc_new/"
branch_dir="/g3/wanghao/kezuo/xhj/GRAPES_GFS3.2/GRAPES_GFS3.2_fix_autobc_allsky/"

remote_dir="${branch_dir}/4DVAR/output/check"
local_dir="/mnt/d/GRAPES/single_obs/${dataTag}/"
checkinno_dir="${local_dir}/checkinno/"
biasprep_dir="${local_dir}/biasprep/"
binary_dir="${local_dir}/binary/"
cdir=$(pwd)

# login and tar data
ssh -i ~/.ssh/id_rsa wanghao@10.40.140.18 \
	"cd ${remote_dir}; \
	tar -cf allsky.tar  ./innotbinno_fy3dmwi*.dat_allsky; \
	tar -cf scatt.tar   ./innotbinno_fy3dmwi*.dat_scatt; \
	tar -cf direct.tar  ./innotbinno_fy3dmwi*.dat_direct; \
	tar -cf obs_checkinno.tar  ./obs_checkinno*.dat; \
	cd ../; \
	tar -cf biasprep.tar biasprepfy3dmwi*.dat_*"

# create local directory
if [ ! -d ${local_dir} ]; then
	mkdir ${local_dir}
fi

if [ ! -d ${binary_dir} ]; then
	mkdir ${binary_dir}
fi

# download data
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/allsky.tar ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/scatt.tar  ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/direct.tar ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/obs_checkinno.tar  ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../biasprep.tar ${local_dir}

scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/mini_info*.dat ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../../std.error.0000 ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../../std.out.0000 ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../../namelist.4dvar ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../../namelist.input ${local_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../../namelist_l.input ${local_dir}

scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../dxa.grd ${binary_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../dxa.ctl  ${binary_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../xa.ctl ${binary_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../xa.grd ${binary_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../xb.ctl ${binary_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../xb.grd ${binary_dir}
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../grapes_input* ${binary_dir}/grapes_input.out
scp -i  ~/.ssh/id_rsa -r wanghao@10.40.140.18:${remote_dir}/../../grapes_input ${binary_dir}/grapes_input.in

# Generate ctl files for grapes_input
cd ${binary_dir} 
rm ./grapes_input.out.ctl ./grapes_input.in.ctl
echo "dset ^grapes_input.out" >> ./grapes_input.out.ctl
echo "dset ^grapes_input.in"  >> ./grapes_input.in.ctl
grapesinputctlBody=$"
undef 2.56E-9\n\
title xa field\n\
FILEHEADER 404\n\
xdef 1440 linear 0 0.25\n\
ydef 721 linear -90 0.25\n\
zdef 89 linear 1 1\n\
tdef 1 linear 00z01NOV2020 6hr\n\
vars 17\n\
pi 89 99 height in m\n\
u 89 99 height in m\n\
v 89 99 height in m\n\
w 89 99 height in m\n\
th 89 99 height in m\n\
qv 89 99 height in m\n\
qc 89 99 height in m\n\
qr 89 99 height in m\n\
qi 89 99 height in m\n\
qs 89 99 height in m\n\
qg 89 99 height in m\n\
pcf 89 99 height in m\n\
ts 0 99 height in m\n\
t2 0 99 height in m\n\
q2 0 99 height in m\n\
u10 0 99 height in m\n\
v10 0 99 height in m\n\
endvars"

echo -e $grapesinputctlBody >> ./grapes_input.out.ctl
echo -e $grapesinputctlBody >> ./grapes_input.in.ctl

# # Get netCDF4 File
cd $cdir
source /home/xhj/software/miniconda3/etc/profile.d/conda.sh
conda activate rdop
python test_convert_to_nc.py "${binary_dir}/grapes_input.out.ctl" grapes_input
python test_convert_to_nc.py "${binary_dir}/grapes_input.in.ctl" grapes_input
python test_convert_to_nc.py "${binary_dir}/dxa.ctl" dxa
python test_convert_to_nc.py "${binary_dir}/xa.ctl" xa
python test_convert_to_nc.py "${binary_dir}/xb.ctl" xa

cd $binary_dir
rm *.ctl
rm *.grd
rm *.in
rm *.out

# Untar
if [ ! -d ${checkinno_dir} ]; then
	mkdir ${checkinno_dir}
fi

if [ ! -d ${biasprep_dir} ]; then
	mkdir ${biasprep_dir}
fi

cd ${local_dir}
tar -xf allsky.tar -C ${checkinno_dir} && rm allsky.tar
tar -xf scatt.tar -C ${checkinno_dir} && rm scatt.tar
tar -xf direct.tar -C ${checkinno_dir} && rm direct.tar
tar -xf obs_checkinno.tar -C ${checkinno_dir} && rm obs_checkinno.tar
tar -xf biasprep.tar -C ${biasprep_dir} && rm biasprep.tar

# link
cd ${cdir}
ln -sf ${local_dir} .

exit
