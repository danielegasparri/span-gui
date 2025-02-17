#script to install the python packages needed for SPAN
echo -e "This script will check and install the Python modules needed to run SPAN\n"
read -n 1 -s -r -p "If you are ok with that, press any key to continue"
pip3 install -r python_packages.txt
echo -e "###############################################\n"
echo -e "Python modules installed. Now you can use SPAN!\n"
echo -e "In the terminal, just do: python3 __main__.py \n"
echo -e "###############################################\n"
