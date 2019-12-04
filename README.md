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

où * est une fonction blas ou lapack  (e.g. axpy, ger, scal, gemm, ...)
(dans tst/ il n'y a que dgemm, dgetrf plus d'autre fonction utiles pour les mesures de performances)

## Environnement GUIX utilisé

** Guix pas disponible pour le moment, vous devez utiliser la commande "module"  comme mentionner plus bas **
Pour avoir le même environnement guix que pour les tests développer :
```sh
guix environment --pure mkl git --ad-hoc less gcc-toolchain coreutils mkl emacs -- /bin/bash
```

En cas d'OPENMP :
```sh
export OMP_NUM_THREADS=40
```

**Warning** : Pour compiler les tests getrf, il ne faut pas se mettre sous *guix*, mais être connecter à Plafrim et faire:
```sh
module load compiler/gcc/9.1.0 compiler/intel/2019_update4
```
Il faut également retirer du *makefile* la version de *MKL* utilisée (conflit avec celle utilisée pour compiler les testeurs de la *libalgonum.so*.

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

### Dgemm
Dgemm-perf  prend 4 arguments en entrée
```sh
./tst/dgemm-perf <start> <end> <step> <nsample>
```
* start   : taille de début
* end     : taille maximale (pas atteinte si il n'y a pas de n tel que end=start*((1 + step/100)^n))
* step    : augmentation de la taille du problème (en %). Ex : si step = 25, alors la les valeurs testées seront de la forme start*(1.25^n)
* nsample : nombre d'échantillon à tester pour chaque taille de matrice traitée.

### Dgemm tuilé

## Nettoyer le projet

Pour nettoyer le projet, faire un 
```sh
make clean
```
