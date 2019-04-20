#include <stdio.h>
#include <math.h>

#define max 20

typedef float Matrix[max][max];
typedef float Vet[max];

void limpa (Matrix A){
     for (int i = 0 ; i<max ; i++){
            for (int j = 0 ; j<max ; j++){
                A[i][j] = 0;
            }
     }
}

void Impm (int n, Matrix A) {

     for (int i = 0 ; i<n ; i++){
            printf ("\n");
            for (int j = 0 ; j<max ; j++){
                printf ("%.3f ", A[i][j]);
            }
        }
}

int convDD (int n, Matrix A){
     float soma;
     for (int i = 0 ; i < n ; i++){
         soma = 0;
         for (int j = 0 ; j < n ; j++){
             if (i!=j){
                soma = soma + fabs(A[i][j]);
             }
         }
         if (A[i][i] < soma){
            return (0);
         }
     }
     return (1);
}

int convVL(int n, Matrix A){
    float soma;
    for (int i=0 ; i < n ; i++){
         soma=0;
         for (int j=0 ; j < n ; j++){
             if (i!=j)
                soma = soma + fabs(A[i][j]/A[i][i]);
         }
         if (soma>1){
            return (0);
         }
     }
     return (1);
}

int convVC(int n, Matrix A){
    float soma;
    for (int i=0 ; i < n ; i++){
         soma=0;
         for (int j=0 ; j < n ; j++){
             if (i!=j)
                soma = soma + fabs(A[j][i]/A[j][j]);
         }
         if (soma>1){
            return (0);
         }
     }
     return (1);
}

int convSas(int n, Matrix A){
     float B[n];
     int VS;
     float soma=0;
     for (int i=0 ; i < n ; i++){
         B[i]=1;
     }

     for (int i=0 ; i < n ; i++){

         if (soma<1){
            soma=0;
            VS=1;
            for (int j=0 ; j < n ; j++)
                if (i!=j)
                   soma = soma + (fabs(A[i][j]/A[i][i]))*B[j];
            if (soma>=1){
               B[i]=soma;
               VS=0;
            }
            B[i]=soma;
         }
     }
     return VS;
}


void SistTriInf(int n,Matrix A,Vet b,Vet y) {
     float soma;
     for (int i = 0 ; i < n ; i++) {
         soma = 0;
         for (int j = 0 ; j < i ; j++)
             soma = soma + y[j]*A[i][j];
         y[i]=(b[i]-soma)/A[i][i];
     }
}

void SistTriSup(int n,Matrix A,Vet y,Vet x) {
     float soma;
     for (int i = n-1 ; i >= 0;i--) {
         soma=0;
         for (int j=i+1 ; j < n ; j++)
             soma = soma + x[j]*A[i][j];
         x[i]=(y[i] - soma)/A[i][i];
     }
}

void DecomposicaoLU(int n,Matrix A,Vet b,Vet x){

     Matrix U,L;
     limpa (U);
     limpa (L);
     Vet y;

     for (int i=0 ; i < n ; i++){
         for (int j=i ; j < n ; j++){
             U[i][j]=A[i][j];
             L[i][i]=1;
             for (int k=0 ; k < i; k++) {
                 U[i][j]=U[i][j] - (L[i][k]*U[k][j]);
             }
         }

         for (int m=i+1 ; m<n ; m++){
             L[m][i] = A[m][i];
             for (int k=0 ; k < i ; k++) {
                 L[m][i]=A[m][i] - (L[m][k]*U[k][i]);
             }
             L[m][i]=L[m][i]/U[i][i];
         }
     }


     printf ("A matriz L eh: \n");
     Impm (max, L);

     printf ("\n\n");
     printf ("A matriz U eh: \n");
     Impm (max, U);

     printf ("\n\n");


     SistTriInf(n,L,b,y);
     SistTriSup(n,U,y,x);

}

int JacobiR(int n , Matrix A , Vet b , Vet x0 , Vet x , int *it) {
         int i,j;
         float eps = 0.00001;
         float dx, maior, erro;
         float soma=0;


         if ((convVL(n, A) || convVC(n, A) || convDD(n, A)) == 0){
            printf ("Nao converge por Jacobi - Richardson");
            return (0);
         }

         for (i=0 ; i < n ; i++) {
             x[i] = x0[i];
         }

         do {
            for (i=0 ; i < n ; i++) {
                x0[i]=x[i];
            }

            for (i=0 ; i < n ; i++){
                x[i] = b[i];
                for (j=0 ; j < n ; j++){
                    if (i!=j)
                       x[i] = x[i] - A[i][j] * x0[j];
                }
                x[i] = x[i]/A[i][i];
            }

            (*it)++;

            dx = 0;
            maior = 0;
            for (i=0 ; i < n ; i++){
                if (fabs(x0[i] - x[i]) > dx){
                   dx = fabs(x[i] - x0[i]);
                }
                if (fabs(x[i]) > maior){
                   maior = x[i];
                }
            }
            erro = dx/maior;
     } while (erro > eps);
}




int GaussSeidel(int n , Matrix A , Vet b , Vet x0 , Vet x , int *it) {
     float erro, dx;
     float eps = 0.00001;
     float soma, maior;

     if ((convSas(n,A) || convVL(n,A) || (convVC(n, A)) || convDD(n, A))==0) {
        printf ("Nao converge por Gauss Seidel");
        return (0);
     }

     for (int i=0 ; i < n ; i++) {
         x[i]=x0[i];
     }

     do {
            for (int i=0 ; i < n ; i++) {
                x0[i]=x[i];
            }

            for (int i=0 ; i < n ; i++){
                x[i]=b[i];
                for (int j=0 ; j < n ; j++)
                    if (i!=j)
                       x[i] = x[i] - A[i][j] * x[j];
                x[i] = x[i]/A[i][i];
            }

            dx=0;
            maior = 0;
            (*it)++;
            for (int i=0 ; i < n ; i++){
                if (fabs(x0[i] - x[i]) > dx) {
                   dx = fabs(x0[i] - x[i]);
                }
                if (fabs(x[i]) > maior) {
                   maior = x[i];
                }
            }
            erro = dx/maior;

     } while (erro > eps);
}

void ImpSol (int n,Vet x){
       printf ("O vetor solucao eh : \n");
        printf ("(");
        for (int j = 0 ; j<max ; j++){
            printf ("%.3f ", x[j]);
        }
        printf (")(t)");
}

int leitura (int n, Matrix A, Vet b){
     int i=0, j=-1;
     float aux;
     FILE *arq;
     if (!(arq = fopen ("matriz.txt" , "r"))){
        return (0);
     }
     while (i<n){
           j++;
           fscanf (arq , "%f" , &aux);
           A[i][j] = aux;
           if (j == 24){
              j = -1;
              i++;
           }
     }
     fclose (arq);
     i = 0;
     arq = fopen ("b.txt" , "r");
     while (i<n){
           fscanf (arq , "%f" , &aux);
           b[i] = aux;
           i++;
     }
     fclose (arq);
     return (1);
}

int main (){

    Matrix A;
    int it = 0;
    Vet b , x , x0;

    for (int i = 0 ; i < max ; i++){
        x0[i] = 0;
    }
    if (leitura (max, A, b)){
        int opc;
        printf ("1 - Decomposicao LU\n");
        printf ("2 - Jacobi\n");
        printf ("3 - Gauss\n >> ");
        scanf ("%d", &opc);
        switch (opc) {
               case 1: DecomposicaoLU(max , A , b , x);
                       break;
               case 2: JacobiR(max , A , b , x0, x, &it);
                       printf ("\nO numero de iteracoes gastas pelo metodo foi: %d\n", it);
                       break;
               case 3: GaussSeidel( max , A , b , x0 , x , &it);
                       printf ("\nO numero de iteracoes gastas pelo metodo foi: %d\n", it);
                       break;
        }

        ImpSol (max,x);
    }
    else {
         printf ("ERRO: O arquivo não pode ser aberto");
    }
}
