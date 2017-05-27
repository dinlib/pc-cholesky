# Cholesky Parallel vs Sequential
Primeiro trabalho da disciplina de programação concorrente
Aluno: Paulo Alexandre Piornedo Panucci   RA: 88380

# Compilando
Para compilar os códigos paralelos e o sequenciais, vá até a pasta src/ no terminal e digite make compile.

#Rodando
Sequencial Polybench:
  ./cholesky.out
Sequencial Desenvolvido:
  ./cholesky_developed.out 1 para rodar com as informações de cache miss do papi
  ./cholesky_developed.out 2 para rodar com as informações de ciclos e instruções do papi
OpenMP:
  ./cholesky_omp.out $numthreads 1 para rodar com as informações de cache miss do papi, e $numthreads o número de threads
  ./cholesky_omp.out $numthreads 2 para rodar com as informações de ciclos e instruções do papi, e $numthreads o número de threads
Pthread:
  ./cholesky_pthread.out $numthreads 1 para rodar com as informações de cache miss do papi, e $numthreads o número de threads
  ./cholesky_pthread.out $numthreads 2 para rodar com as informações de ciclos e instruções do papi, e $numthreads o número de threads

Para remover os executáveis: make clean.

#Scripts
  o script run.sh foi utilizado para rodar a análise de speedup, e para utilizálo é importate tirar todos os prints presentes no código, deixando somente o print de tempo de execução. Além disso para o seu funcionamento deve-se ter uma pasta data/ no mesmo local do script e dos executáveis.

  o script papi.sh foi utilizado para rodar a análise do papi, e paara o seu funcionamento deve-se ter uma pasta papidata/ no mesmo local do script e dos executáveis.
