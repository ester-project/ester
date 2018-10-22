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

        stage('Test GNU') {
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

        stage('Test Clang') {
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
    }

    post {
        always {
            dir('build-gnu') {
                deleteDir()
            }
            dir('build-clang') {
                deleteDir()
            }
        }
        failure {
            mail to: 'bertrand.putigny@gmail.com',
                 subject: "Failed Jenkins pipeline: ${currentBuild.fullDisplayName}",
                 body: "Something is wrong with ${env.BUILD_URL}"
        }
    }
}
