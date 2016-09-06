/home/marpodbuser/marpodbENV/etc/marodb-express/apachectl stop
rm -rf /home/marpodbuser/marpodbENV/etc/marodb-express
mod_wsgi-express setup-server marpodb.wgsi --port=8081 --server-root=/home/marpodbuser/marpodbENV/etc/marodb-express
/home/marpodbuser/marpodbENV/etc/marodb-express/apachectl start
