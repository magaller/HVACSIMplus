copy hvacsim.ctf thrizone.ctf
copy hvacsim.met thrizone.met
modsim < path_1_thrizone.inp
copy thrizone.fin thrizone.ini
modsim < path_2_thrizone.inp
sortsb < path_3_thrizone.inp

