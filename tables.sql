DROP SCHEMA public CASCADE;
CREATE SCHEMA public;

CREATE TABLE gene(
	id serial,
	name varchar(100),
	seq text,
	PRIMARY KEY(id)
);

CREATE TABLE transcript(
	id serial,
	name varchar(100),
	gene_id int REFERENCES gene,
	coordinates text,
	PRIMARY KEY(id)
); 

CREATE TABLE cds(
	id serial,
	name varchar(100),
	transcript_id int REFERENCES transcript,
	seq text,
	coordinates text,
	PRIMARY KEY(id)
);

CREATE TABLE pfam_hit(
	id serial,
	cds_id int REFERENCES cds,
	name varchar(100),
	acc varchar(100),
	e_val double precision,
	c_val double precision,
	descr text,
	coordinates text,
	PRIMARY KEY(id)
);

CREATE TABLE blastp_hit(
	id serial,
	cds_id int REFERENCES cds,
	uni_id varchar(100),
	protein_name text,
	gene_name text,
	origin text,
	e_val double precision,
	coordinates text,
	qlen int,
	tlen int
);

CREATE TABLE blastm_hit(
	id serial,
	transcript_id int REFERENCES transcript,
	gb_id varchar(100),
	protein_name text,
	gene_name text,
	origin text,
	e_val double precision,
	coordinates text,
	qlen int,
	tlen int
);
