#!groovy

pipeline {
    agent any

        triggers {
            pollSCM('H/15 * * * *')
        }

    stages {

        stage('Bootstrap') {
            steps {
                sh './bootstrap'
            }
        }

        stage('Autotools build with GNU') {
            steps {
                sh 'mkdir -p build-gnu'
                dir('build-gnu') {
                    sh "../configure CC=gcc CXX=g++ FC=gfortran --prefix=\$(pwd)/"
                    sh "make -j4 && make install"
                    sh "./bin/ester 1d -noplot -o M5-1d.h5"
                    sh "./bin/ester 2d -noplot -i M5-1d.h5 -o M5-omega0.8.h5 -Omega_bk 0.8"
                }
            }
        }
        stage('CMake build with GNU') {
            steps {
                sh 'mkdir -p cmake-build-gnu'
                dir('cmake-build-gnu') {
                    sh "cmake -DCMAKE_INSTALL_PREFIX=\$PWD .."
                    sh "make -j4 && make install"
                    sh "./bin/ester 1d -noplot -o M5-1d.h5"
                    sh "./bin/ester 2d -noplot -i M5-1d.h5 -o M5-omega0.8.h5 -Omega_bk 0.8"
                }
            }
        }

        stage('Autotools build with Clang') {
            steps {
                sh 'mkdir -p build-clang'
                dir('build-clang') {
                    sh "../configure CC=clang CXX=clang++ FC=gfortran --prefix=\$(pwd)/"
                    sh "make -j4 && make install"
                    sh "./bin/ester 1d -noplot -o M5-1d.h5"
                    sh "./bin/ester 2d -noplot -i M5-1d.h5 -o M5-omega0.8.h5 -Omega_bk 0.8"
                }
            }
        }
        stage('CMake build with Clang') {
            steps {
                sh 'mkdir -p cmake-build-clang'
                dir('build-clang') {
                    sh "cmake -DCMAKE_INSTALL_PREFIX=\$PWD -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang .."
                    sh "make -j4 && make install"
                    sh "./bin/ester 1d -noplot -o M5-1d.h5"
                    sh "./bin/ester 2d -noplot -i M5-1d.h5 -o M5-omega0.8.h5 -Omega_bk 0.8"
                }
            }
        }

        stage('Read reference') {
            steps {
                dir('build-gnu') {
                    sh "./bin/ester info ../references/M5-omega0.8.h5"
                    sh "./bin/ester vtk ../references/M5-omega0.8.h5 -o M5-omega0.8.vtk"
                }
            }
        }

    }

    post {
        always {
            dir('build-gnu') {
                deleteDir()
            }
            dir('build-clang') {
                deleteDir()
            }
            dir('cmake-build-clang') {
                deleteDir()
            }
            dir('cmake-build-clang') {
                deleteDir()
            }
        }
    }
}
