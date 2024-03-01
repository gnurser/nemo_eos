cd ~/Dropbox/VC_code/nemo_eos
mamba activate big
site_packages=$(get_site_packages) # => /Users/agn/miniforge3/envs/big/python3.10/site-packages
mamba env create -f  nemo_eos_dev.yml
mamba activate nemo_eos_dev
meson setup builddir -Dprefix=/Users/agn/miniforge3/envs/big
meson compile -C builddir
meson install -C builddir
# Note *do not use* rest of $site_packages i.e python3.10/site-packages bit
