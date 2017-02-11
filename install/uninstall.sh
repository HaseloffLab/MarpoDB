#!/bin/bash
echo    ""
read -p "## Are you sure you want to uninstall?" -n 1 -r
echo    ""
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
fi

cd ..
rm ~/.paths
rm ~/.pypaths
rm ~/.ldpaths
rm -Rf ~/var/
rm -Rf  ~/ENV/
cd ..
rm -Rf MarpoDB/