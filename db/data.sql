  begin;

  delete from user_account;

  copy user_account (id,email,email_confirmed,created_at,updated_at) from stdin (format csv);
"61b64c10-b3bc-11e7-9b1c-9f971ce2bd8c","zhuxun2@gmail.com","true","2017-12-04 16:04:10.969000","2017-12-04 16:04:10.969000"
"e262c322-df30-4f21-8ecb-464a1b940230","opoirion@cc.hawaii.edu","false","2017-12-04 16:04:10.972000","2017-12-04 16:04:10.972000"
"2c4045e7-f806-42a8-970c-363c8b41c3dd","falakwaa@cc.hawaii.edu","false","2017-12-04 16:04:10.973000","2017-12-04 16:04:10.973000"
"017617ff-fd29-4ed4-b6c8-8f9a25bca9a1","pbenny@cc.hawaii.edu","false","2017-12-04 16:04:10.973000","2017-12-04 16:04:10.973000"
"25afb4ca-6534-4d41-aac5-01a9b9ee1ca1","CLassiter@cc.hawaii.edu","false","2017-12-04 16:04:10.973000","2017-12-04 16:04:10.973000"
"0205bbb9-b533-49c1-a9ff-a3244747eac1","byunits@cc.hawaii.edu","false","2017-12-04 16:04:10.973000","2017-12-04 16:04:10.973000"
"e74dae06-8ac4-4bd8-a346-13cb37e44c51","twolfgruber@cc.hawaii.edu","false","2017-12-04 16:04:10.973000","2017-12-04 16:04:10.973000"
"c7f4e710-20d6-4837-b194-6517b9e797ca","carisdak@hawaii.edu","false","2017-12-04 16:04:10.973000","2017-12-04 16:04:10.973000"
\.

  delete from user_profile;

  copy user_profile (user_id,display_name,gender,institution,picture,created_at,updated_at) from stdin (format csv);
61b64c10-b3bc-11e7-9b1c-9f971ce2bd8c,Xun Zhu,male,"University of Hawaii, Manoa",http://garmiregroup.org/member_photos/xun_zhu.jpg,2017-11-30 15:20:09.943000,2017-11-30 15:20:09.943000
e262c322-df30-4f21-8ecb-464a1b940230,Olivier Poirion,male,University of Hawaii Cancer Center,http://garmiregroup.org/member_photos/olivier_poirion.jpg,2017-11-30 15:20:09.948000,2017-11-30 15:20:09.948000
2c4045e7-f806-42a8-970c-363c8b41c3dd,Fadhl Alakwaa,male,University of Hawaii Cancer Center,http://garmiregroup.org/member_photos/fadhl_alakwaa.jpg,2017-11-30 15:20:09.954000,2017-11-30 15:20:09.954000
017617ff-fd29-4ed4-b6c8-8f9a25bca9a1,Paula Benny,female,University of Hawaii Cancer Center,http://garmiregroup.org/member_photos/paula_benny.jpg,2017-11-30 15:20:09.954000,2017-11-30 15:20:09.954000
25afb4ca-6534-4d41-aac5-01a9b9ee1ca1,Cameron Lassiter,male,University of Hawaii Cancer Center,http://garmiregroup.org/member_photos/cameron_lassiter.png,2017-11-30 15:20:09.954000,2017-11-30 15:20:09.954000
0205bbb9-b533-49c1-a9ff-a3244747eac1,Breck Yunits,male,University of Hawawii Cancer Center,http://garmiregroup.org/member_photos/breck_yunits.png,2017-11-30 15:20:09.954000,2017-11-30 15:20:09.954000
e74dae06-8ac4-4bd8-a346-13cb37e44c51,Thomas Wolfgruber,male,University of Hawaii Cancer Center,http://garmiregroup.org/member_photos/thomas_wolfgruber.jpg,2017-11-30 15:20:09.954000,2017-11-30 15:20:09.954000
c7f4e710-20d6-4837-b194-6517b9e797ca,Cedric Arisdakessian,male,University of Hawaii Cancer Center,http://garmiregroup.org/member_photos/cedric_arisdakessian.jpg,2017-11-30 15:20:09.954000,2017-11-30 15:20:09.954000
\.

  delete from project;

  copy project (id,owner_id,rank,name,description,created_at,updated_at) from stdin (format csv);
\.

  delete from step;

  copy step (id,project_id,status,rank,gbox,args,created_at,updated_at) from stdin (format csv);
\.

  end;
