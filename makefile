all:	SETL genWorld SETL_par

SETL:	SETL.c
	gcc -o SETL SETL.c

genWorld:	genWorld.c
	gcc -o genWorld genWorld.c

SETL_par: SETL_par.c
	mpicc -o SETL_par SETL_par.c
