if [ "$#" -eq 0 ]; then
	wget -N https://www.deryase.org/Downloads/Installer/deryaSE_installer_rel_58.bin
elif [ "$1" == "--update" ]; then
	\rm -f deryaSE_installer_rel_58.bin
	wget https://www.deryase.org/Downloads/Installer/deryaSE_installer_rel_58.bin
else
	echo "Unknown parameter [$1]"
	exit 1
fi

chmod 0755 deryaSE_installer_rel_58.bin
./deryaSE_installer_rel_58.bin --acceptLicence --local

