# /usr/bin/env ruby

begin
	a = File.read(ARGV.shift)
	b = a.split('<div class="grid_19 omega main-content">')
	c = b[1].split("</body>")
	d = ' <div class="container_24"><div class="grid_24 clearfix" id="content"><div class="grid_19 omega main-content">'+c[0]
	e = d.gsub('resources/' , '../static/interpro/resources/')
	f = e.gsub('pfam.sanger.ac.uk' , 'pfam.xfam.org')
	puts f
end
