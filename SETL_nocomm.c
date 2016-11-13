#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
/*
MPI Global Variables
*/
int slaves;
int myid;
//#define DEBUG
#define MASTER_ID slaves
/***********************************************************
  Helper functions 
***********************************************************/

//For exiting on error condition
void die(int lineNo);

//For trackinng execution
long long wallClockTime();


/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue );

char** allocateSquareMatrixNoEmpty( int size, char* tmpChars );

char** allocateMatrixNoEmpty( int row, int col, char* tmpChars );

char** allocateMatrix( int row, int col, char defaultValue );

void freeSquareMatrix( char** );

void printSquareMatrix( char**, int size );


/***********************************************************
   World  related functions
***********************************************************/

#define ALIVE 'X' 
#define DEAD 'O'

char** readWorldFromFile( char* fname, int* size );

int countNeighbours(char** world, int row, int col);

void evolveWorld(char** curWorld, char** nextWorld, int row, int col);


/***********************************************************
   Simple circular linked list for match records
***********************************************************/

typedef struct MSTRUCT {
    int iteration, row, col, rotation;
    struct MSTRUCT *next;
} MATCH;


typedef struct {
    int nItem;
    MATCH* tail;
} MATCHLIST;

MATCHLIST* newList();

void deleteList( MATCHLIST*);

void insertEnd(MATCHLIST*, int, int, int, int);

void printList(MATCHLIST*);

int matchToInt(MATCH *mat);

int* transferListToArr(MATCHLIST* list);
MATCH* intToMatch(int index, int iteration);
/***********************************************************
   Search related functions
***********************************************************/

//Using the compass direction to indicate the rotation of pattern
#define N 0 //no rotation
#define E 1 //90 degree clockwise
#define S 2 //180 degree clockwise
#define W 3 //90 degree anti-clockwise
#define MAX_FIND_ONCE 100000
char** readPatternFromFile( char* fname, int* size );

void rotate90(char** current, char** rotated, int size);

void searchPatterns(char** world, int wRow, int wCol, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list, int rowOffset);

void searchSinglePattern(char** world, int wSizeRow, int wSizeCol, int interation,
        char** pattern, int pSize, int rotation, MATCHLIST* list, int rowOffset);

int min(int a, int b){
    if (a < b) return a; else return b;
}

int sortFunction( const void *a, const void *b);
/***********************************************************
   Main function
***********************************************************/

int masterWork(int argc, char** argv){
    char **curW, **nextW, **temp, dummy[20];
    char **patterns[4];
    int dir, iterations, iter;
    int size, patternSize;
    long long before, after;
    MATCHLIST* list, *tmpList;
    MPI_Status Stat;
    int sendTag = 0;
    if (argc < 4 ){
        fprintf(stderr, 
            "Usage: %s <world file> <Iterations> <pattern file>\n", argv[0]);
        exit(1);
    } 

    curW = readWorldFromFile(argv[1], &size);
    nextW = allocateSquareMatrix(size+2, DEAD);

    //Start timer
    before = wallClockTime();

    printf("World Size = %d\n", size);

    iterations = atoi(argv[2]);
    printf("Iterations = %d\n", iterations);

    patterns[N] = readPatternFromFile(argv[3], &patternSize);
    for (dir = E; dir <= W; dir++){
        patterns[dir] = allocateSquareMatrix(patternSize, DEAD);
        rotate90(patterns[dir-1], patterns[dir], patternSize);
    }
    printf("Pattern size = %d\n", patternSize);

    /*Send size and iteration information all slaves*/
    int basicInfo[3] = {size, iterations, patternSize};
    for (int i = 0; i < slaves; i++){
        MPI_Send(basicInfo, 3, MPI_INT, i, sendTag, MPI_COMM_WORLD);
    }


#ifdef DEBUG
    printSquareMatrix(patterns[N], patternSize);
    printSquareMatrix(patterns[E], patternSize);
    printSquareMatrix(patterns[S], patternSize);
    printSquareMatrix(patterns[W], patternSize);
#endif
    for (int i = N; i <= W; i++ ){
        sendTag++;
        for (int j = 0; j < slaves; j++){
            MPI_Send(patterns[i][0], patternSize * patternSize, MPI_CHAR, j, sendTag, MPI_COMM_WORLD);
        }
    }    
    
    sendTag++;
    int responsibleRows[slaves];
    for (int i = 0; i < slaves; i++){
        responsibleRows[i] = size / slaves;
        if (i < size % slaves) responsibleRows[i]++;
    }
    int currentRow = 1; //Start from row 1 as row 0 is meaningless
    for (int i = 0; i < slaves; i++){
        int stopRow = min(currentRow + responsibleRows[i] -1 + patternSize-1, size+1); //stops at size row as this is the last meaningful row
        int tmpSize = (stopRow - currentRow +2) * (size+2);
        MPI_Send(curW[currentRow-1], tmpSize, MPI_CHAR, i, sendTag, MPI_COMM_WORLD);
        currentRow += responsibleRows[i];
    }



    //Actual work start
    list = newList();
    for (iter = 0; iter < iterations; iter++){
        tmpList = newList();
        for (int i = 0; i < slaves; i++){
            int matchSize;
            MPI_Recv(&matchSize, 1, MPI_INT, i, iter, MPI_COMM_WORLD, &Stat);
            int* matchArr = (int *) malloc(sizeof(int) * matchSize);
            MPI_Recv(matchArr, matchSize, MPI_INT, i, iter, MPI_COMM_WORLD, &Stat);
            for (int j = 0; j < matchSize; j++){
                MATCH *newMatch = intToMatch(matchArr[j], iter);
                insertEnd(tmpList, newMatch->iteration, newMatch->row, newMatch->col, newMatch->rotation);
            }
        }
        int* tmpArr = transferListToArr(tmpList);

        qsort((void *)tmpArr, tmpList->nItem, sizeof(int), sortFunction);
        for (int j = 0; j < tmpList->nItem; j++){
            MATCH *newMatch = intToMatch(tmpArr[j], iter);
            insertEnd(list, newMatch->iteration, newMatch->row, newMatch->col, newMatch->rotation);
            
        }            
    }
//     for (iter = 0; iter < iterations; iter++){

// #ifdef DEBUG
//         printf("World Iteration.%d\n", iter);
//         printSquareMatrix(curW, size+2);
// #endif

//         searchPatterns( curW, size, iter, patterns, patternSize, list);

//         //Generate next generation
//         evolveWorld( curW, nextW, size );
//         temp = curW;
//         curW = nextW;
//         nextW = temp;
//     }


    printList( list );

    //Stop timer
    after = wallClockTime();

    printf("Parallel SETL took %1.2f seconds\n", 
        ((float)(after - before))/1000000000);


//     //Clean up
//     deleteList( list );

//     freeSquareMatrix( curW );
//     freeSquareMatrix( nextW );

//     freeSquareMatrix( patterns[0] );
//     freeSquareMatrix( patterns[1] );
//     freeSquareMatrix( patterns[2] );
//     freeSquareMatrix( patterns[3] );

}

int slaveWork(){
    char **patterns[4];
    int basicInfo[3];
    int size, patternSize, iterations;
    int receiveTag = 0;
    char **curW, **nextW, **temp;
    MPI_Status status;
    MATCHLIST* list;



    list = newList();
    MPI_Recv(basicInfo, 3, MPI_INT, MASTER_ID, receiveTag, MPI_COMM_WORLD, &status);
    size = basicInfo[0];
    iterations = basicInfo[1];
    patternSize = basicInfo[2];
#ifdef DEBUG
    printf("Slave node %d received size = %d iterations = %d patternSize = %d\n", myid, size, iterations, patternSize);
#endif
    char tmpChars[4][patternSize * patternSize];
    for (int i = 0; i < 4; i++){
        receiveTag++;
        //patterns[dir] = allocateSquareMatrix(patternSize, DEAD);
        MPI_Recv(tmpChars[i], patternSize * patternSize, MPI_CHAR, MASTER_ID, receiveTag, MPI_COMM_WORLD, &status);
        patterns[i] = allocateSquareMatrixNoEmpty(patternSize, tmpChars[i]);
    }
#ifdef DEBUG
    printf("Slave node %d received four pattern matrix\n", myid);
    printSquareMatrix(patterns[N], patternSize);
    printSquareMatrix(patterns[E], patternSize);
    printSquareMatrix(patterns[S], patternSize);
    printSquareMatrix(patterns[W], patternSize);
#endif
    
    receiveTag++;
    int responsibleRows[slaves];
    for (int i = 0; i < slaves; i++){
        responsibleRows[i] = size / slaves;
        if (i < size % slaves) responsibleRows[i]++;
    }
    int currentRow = 1; //Start from row 1 as row 0 is meaningless
    char matrixInfo[(responsibleRows[myid] + patternSize) * (size+2)];
    int myRowNumber;
    int rowOffset;
    for (int i = 0; i < slaves; i++){
        int stopRow = min(currentRow + responsibleRows[i] -1 + patternSize-1, size+1); //stops at size row as this is the last meaningful row
        int tmpSize = (stopRow - currentRow +2) * (size+2);
        if (i == myid){
            MPI_Recv(matrixInfo, tmpSize, MPI_CHAR, MASTER_ID, receiveTag, MPI_COMM_WORLD, &status);
            myRowNumber = stopRow - currentRow +2;
            rowOffset = currentRow - 1;
            break;
        }
        currentRow += responsibleRows[i];
    }
    curW = allocateMatrixNoEmpty((size + 2), myRowNumber, matrixInfo);
    nextW = allocateMatrix((size + 2), myRowNumber, DEAD);
#ifdef DEBUG
    for (int i = 1; i < myRowNumber; i++){
        for (int j = 1; j <= size; j++){
            printf("%c",curW[i][j]);
        }
        printf("\n");
    }
#endif
    //searchPatterns( curW, myRowNumber-1, size, 0, patterns, patternSize, list, rowOffset);
    //printList(list);
    int sendTag = 0;
    char buffer[size];
    for (int i = 0; i< iterations; i++){
        searchPatterns( curW, myRowNumber-1, size, i, patterns, patternSize, list, rowOffset);
        evolveWorld(curW, nextW, myRowNumber-2, size);
        temp = curW;
        curW = nextW;
        nextW = temp;
        
#ifdef DEBUG
        if (myid == 1 && i == 1){
            printf("world is like!\n");
            for (int q = 0; q < myRowNumber; q++){
                for (int p = 1; p <= size; p++){
                    printf("%c", curW[q][p]);
                }
                printf("\n");
            }
        } 
#endif
        // if (myid != 0){
        //     for (int j = 0; j < patternSize-1; j++){
        //         for (int k = 1; k <= size; k++){
        //             buffer[k-1] = curW[j+1][k];
        //         }
        //         MPI_Send(buffer, size, MPI_CHAR, myid - 1, i * size + myid +j, MPI_COMM_WORLD);
        //     }

        //}
        if (myid != slaves-1){
            for (int k = 1; k <= size/2; k++){
                buffer[k-1] = curW[myRowNumber - patternSize][k];
            }
            MPI_Send(buffer, size/2, MPI_CHAR, myid + 1, i * size +myid, MPI_COMM_WORLD);
        }
        if (myid != 0){
            MPI_Recv(buffer, size/2, MPI_CHAR, myid - 1, i * size + (myid-1), MPI_COMM_WORLD,&status);

            for (int k = 1; k <= size/2; k++){
                curW[0][k] = buffer[k-1];
            }
        }

        // if (myid != slaves-1){
        //     for (int j = 0; j < patternSize-1; j++){
        //         MPI_Recv(buffer, size, MPI_CHAR, myid + 1, i * size + (myid+1) + j, MPI_COMM_WORLD, &status);
        //         for (int k = 1; k <= size; k++){
        //             curW[j + myRowNumber - patternSize +1][k] = buffer[k-1];
        //         }
        //     }            
        // }
        // if (myid == 1 && i == 1){
        //     printf("world after change!\n");
        //     for (int q = 0; q < myRowNumber; q++){
        //         for (int p = 1; p <= size; p++){
        //             printf("%c", curW[q][p]);
        //         }
        //         printf("\n");
        //     }
        // } 

        /*After evolve, transfer the information to neighbours*/
        int *matchArr = transferListToArr(list);
        int matchSize = list->nItem;
        MPI_Send(&matchSize, 1, MPI_INT, MASTER_ID , i, MPI_COMM_WORLD);
        MPI_Send(matchArr, list->nItem, MPI_INT, MASTER_ID , i, MPI_COMM_WORLD);
        list = newList();    
    }
    //printList(list);

}

int main( int argc, char** argv)
{
    int nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    slaves = nprocs - 1;
    if (myid == MASTER_ID){
        masterWork(argc, argv);
    }else{
        slaveWork();
    }
    
    MPI_Finalize();
    return 0;
}

/***********************************************************
  Helper functions 
***********************************************************/


void die(int lineNo)
{
    fprintf(stderr, "Error at line %d. Exiting\n", lineNo);
    exit(1);
}

long long wallClockTime( )
{
#ifdef __linux__
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}

/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue )
{

    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = (char*) malloc(sizeof(char) * size * size);
    if (contiguous == NULL) 
        die(__LINE__);


    memset(contiguous, defaultValue, size * size );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * size );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < size; i++){
        matrix[i] = &contiguous[i*size];
    }

    return matrix;
}
char** allocateMatrix( int row, int col, char defaultValue )
{

    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = (char*) malloc(sizeof(char) * row * col);
    if (contiguous == NULL) 
        die(__LINE__);


    memset(contiguous, defaultValue, col * row );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * col );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < col; i++){
        matrix[i] = &contiguous[i*row];
    }

    return matrix;
}


char** allocateSquareMatrixNoEmpty( int size, char* tmpChars )
{

    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = tmpChars;


    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * size );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < size; i++){
        matrix[i] = &contiguous[i*size];
    }

    return matrix;
}

char** allocateMatrixNoEmpty( int row, int col,  char* tmpChars )
{

    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = tmpChars;


    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * col );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < col; i++){
        matrix[i] = &contiguous[i*row];
    }

    return matrix;
}

void printSquareMatrix( char** matrix, int size )
{
    int i,j;
    
    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            printf("%c", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeSquareMatrix( char** matrix )
{
    if (matrix == NULL) return;

    free( matrix[0] );
}

/***********************************************************
   World  related functions
***********************************************************/

char** readWorldFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **world;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    //Using the "halo" approach
    // allocated additional top + bottom rows
    // and leftmost and rightmost rows to form a boundary
    // to simplify computation of cell along edges
    world = allocateSquareMatrix( size + 2, DEAD );

    for (i = 1; i <= size; i++){
        for (j = 1; j <= size; j++){
            fscanf(inf, "%c", &world[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }

    *sizePtr = size;    //return size
    return world;
    
}

int countNeighbours(char** world, int row, int col)
//Assume 1 <= row, col <= size, no check 
{
    int i, j, count;

    count = 0;
    for(i = row-1; i <= row+1; i++){
        for(j = col-1; j <= col+1; j++){
            count += (world[i][j] == ALIVE );
        }
    }

    //discount the center
    count -= (world[row][col] == ALIVE);

    return count;

}

void evolveWorld(char** curWorld, char** nextWorld, int row ,int col)
{
    int i, j, liveNeighbours;

    for (i = 1; i <= row; i++){
        for (j = 1; j <= col; j++){
            //printf("%d %d\n", i,j);
            liveNeighbours = countNeighbours(curWorld, i, j);
            nextWorld[i][j] = DEAD;
            //Only take care of alive cases
            if (curWorld[i][j] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;
        } 
    }
}

/***********************************************************
   Search related functions
***********************************************************/

char** readPatternFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **pattern;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    pattern = allocateSquareMatrix( size, DEAD );

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            fscanf(inf, "%c", &pattern[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }
    
    *sizePtr = size;    //return size
    return pattern;
}


void rotate90(char** current, char** rotated, int size)
{
    int i, j;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            rotated[j][size-i-1] = current[i][j];
        }
    }
}

void searchPatterns(char** world, int wRow, int wCol, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list, int rowOffset)
{
    int dir;

    for (dir = N; dir <= W; dir++){
        searchSinglePattern(world, wRow, wCol, iteration, 
                patterns[dir], pSize, dir, list, rowOffset);
    }

}

void searchSinglePattern(char** world, int wSizeRow, int wSizeCol, int iteration,
        char** pattern, int pSize, int rotation, MATCHLIST* list, int rowOffset)
{
    int wRow, wCol, pRow, pCol, match;


    for (wRow = 1; wRow <= (wSizeRow-pSize+1); wRow++){
        for (wCol = 1; wCol <= (wSizeCol-pSize+1); wCol++){
            match = 1;
#ifdef DEBUGMORE
            printf("S:(%d, %d)\n", wRow-1, wCol-1);
#endif
            for (pRow = 0; match && pRow < pSize; pRow++){
                for (pCol = 0; match && pCol < pSize; pCol++){
                    if(world[wRow+pRow][wCol+pCol] != pattern[pRow][pCol]){
#ifdef DEBUGMORE
                        printf("\tF:(%d, %d) %c != %c\n", pRow, pCol,
                            world[wRow+pRow][wCol+pCol], pattern[pRow][pCol]);
#endif
                        match = 0;    
                    }
                }
            }
            if (match && (wRow-1 + rowOffset <= wSizeCol-pSize)){

                insertEnd(list, iteration, wRow-1 + rowOffset, wCol-1, rotation);
#ifdef DEBUGMORE
printf("*** Row = %d, Col = %d\n", wRow-1, wCol-1);
#endif
            }
        }
    }
}

/***********************************************************
   Simple circular linked list for match records
***********************************************************/

MATCHLIST* newList()
{
    MATCHLIST* list;

    list = (MATCHLIST*) malloc(sizeof(MATCHLIST));
    if (list == NULL)
        die(__LINE__);

    list->nItem = 0;
    list->tail = NULL;

    return list;
}

void deleteList( MATCHLIST* list)
{
    MATCH *cur, *next;
    int i;
    //delete items first

    if (list->nItem != 0 ){
        cur = list->tail->next;
        next = cur->next;
        for( i = 0; i < list->nItem; i++, cur = next, next = next->next ) {
            free(cur); 
        }

    }
    free( list );
}

void insertEnd(MATCHLIST* list, 
        int iteration, int row, int col, int rotation)
{
    MATCH* newItem;

    newItem = (MATCH*) malloc(sizeof(MATCH));
    if (newItem == NULL)
        die(__LINE__);

    newItem->iteration = iteration;
    newItem->row = row;
    newItem->col = col;
    newItem->rotation = rotation;

    if (list->nItem == 0){
        newItem->next = newItem;
        list->tail = newItem;
    } else {
        newItem->next = list->tail->next;
        list->tail->next = newItem;
        list->tail = newItem;
    }

    (list->nItem)++;

}

void printList(MATCHLIST* list)
{
    int i;
    MATCH* cur;

    printf("List size = %d\n", list->nItem);    


    if (list->nItem == 0) return;

    cur = list->tail->next;
    for( i = 0; i < list->nItem; i++, cur=cur->next){
        printf("%d:%d:%d:%d\n", 
                cur->iteration, cur->row, cur->col, cur->rotation);
    }
}

int* transferListToArr(MATCHLIST* list)
{
    int i;
    MATCH* cur;

    int *arr = (int*) malloc(sizeof(int) * list->nItem);

    if (list->nItem == 0) return arr;

    cur = list->tail->next;
    for( i = 0; i < list->nItem; i++, cur=cur->next){
        //printf("%d:%d:%d:%d\n", cur->iteration, cur->row, cur->col, cur->rotation);
        arr[i] = matchToInt(cur);
    }
    return arr;
}

int matchToInt(MATCH *mat){
    return (mat->row * 10000 + mat->col*1 + mat->rotation*100000000);
} 

MATCH* intToMatch(int index, int iteration){
    MATCH* newItem;
    newItem = (MATCH*) malloc(sizeof(MATCH));
    newItem->col = index % 10000;
    index /= 10000;
    newItem->row = index % 10000;
    index /= 10000;
    newItem->iteration = iteration;
    newItem->rotation = index;
    return newItem;
}

int sortFunction( const void *a, const void *b){
        if(*(int*)a>*(int*)b)
                return 1;
        else if(*(int*)a<*(int*)b)
                return -1;
        else
                return 0;
}