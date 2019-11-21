# Bienvenue dans notre version de LAPACK

Vous trouverez ici les différentes fonctions LAPACK que nous avons implémenté ainsi que des fichier de tests

# Fichiers

 ```
├── Makefile
├── *.o
├── include/
│   	    └── *.h
├── src/
│   	└── *.c
├── tst/
│       ├── *-perf
│   	└── *-perf.c
└── lib/
	└── libmyblas.a
```

where * is a blas or lapack function function (e.g. axpy, ger, scal, gemm, ...)

## Environnement GUIX utilisé

Pour avoir le même environnement guix que pour les tests développer :
```sh
guix environment --pure mkl git --ad-hoc less gcc-toolchain coreutils mkl emacs -- /bin/bash
```

**Warning** : Pour compiler les tests getrf, il ne faut pas se mettre sous *guix*, mais être connecter à Plafrim et faire:
```sh
module load compiler/gcc/9.1.0 compiler/intel/2019_update4
```

## Compiler

Pour compiler le projet faire un 

```sh
make
```

## Test de performance

Pour compiler un fichier en vue de tester ces performance 
```sh
make **nom de la fonction**-perf
```
Puis

### Ddot
Ddot-perf prend 4 arguments en entrée
```sh
./tst/ddot-perf <start> <end> <step> <nsample>
```
* start   : taille de début
* end     : taille maximale (pas atteinte si il n'y a pas de n tel que end=start*((1 + step/100)^n))
* step    : augmentation de la taille du problème (en %). Ex : si step = 25, alors la les valeurs testées seront de la forme start*(1.25^n)
* nsample : nombre d'échantillon à tester pour chaque taille de matrice traitée.

### Dgemm
Dgemm-perf  prend 4 arguments en entrée
```sh
./tst/dgemm-perf <start> <end> <step> <nsample>
```
Ces valeurs ont la même signification que pour ddot.

## Nettoyer le projet

Pour nettoyer le projet, faire un 
```sh
make clean
```
