DROP SCHEMA public CASCADE;
CREATE SCHEMA public;

CREATE TABLE "user"(
	id serial,
	username varchar(50),
	password text,
	email varchar(120),
	is_active boolean,
	affiliation text,
	PRIMARY KEY(id)
);

CREATE TABLE "star_gene"(
	id serial,
	userid int REFERENCES "user",
	geneid int,
	PRIMARY KEY(id)
);
