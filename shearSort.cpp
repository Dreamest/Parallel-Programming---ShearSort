//Studens: Guy & Chen !!!

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <cstdlib>
#include <time.h> //check

//--constants:--
#define ROOT 0
#define ROW_DIR 1
#define COL_DIR 0
#define N_DIMS 2
#define POINT_NUM_OF_ATTRIBUTES 3 // 3 for x,y,z
#define FILE_PATH "C:\\inputFile.txt" //path of input file

//--Point struct:--
typedef struct
{
	double x;
	double y;
	double z;
}Point;

//--------------------------functions decleration:----------------------------//
MPI_Comm createCartComm(int rank, int numOfProcesses, int *coords);

void shearSort(int* coor, Point* myPoint, MPI_Datatype PointType, MPI_Comm cartComm, int pointsInRow);

void loadPointsMatrixFromFile(Point* mat, FILE* fp, int numOfPoints);

double calculateDistanceOfPoint(Point* p1);

void printPointsMatrix(Point* mat, int pointsInRow);

void printAsArray(Point* mat, int pointsInRow);

void createPointType(MPI_Datatype* PointType);


//------------------------------------main:------------------------------------//
int main(int argc, char* argv[])
{
	const int DIMENSIONS = N_DIMS;
	FILE* fp;
	int myId, numOfProcesses;
	MPI_Datatype PointType;
	Point *mat, myPoint;
	int coor[DIMENSIONS];
	int numOfPoints;
	int rowSize;

	MPI_Comm cartComm;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	createPointType(&PointType); //create PointType
	cartComm = createCartComm(myId, numOfProcesses, coor); //initialize cartComm

	//read points matrix from file:
	if (myId == ROOT) //only root read from file
	{
		fopen_s(&fp, FILE_PATH, "r");
		if (!fp) //if file doesn't open - abort.
		{
			printf("File could not open.\n");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}

		fscanf_s(fp, "%d", &numOfPoints); // read numOfPoints from file
		rowSize = (int)sqrt(numOfPoints);

		mat = (Point*)calloc(numOfPoints, sizeof(Point)); //initialize array for points
		if (!mat)
			MPI_Abort(MPI_COMM_WORLD, 2);

		loadPointsMatrixFromFile(mat, fp, numOfPoints); //read all points from file and load them to mat
		fclose(fp); // at this point file is not needed anymore, hence it's being closed.

		//check if num of proccesses isn't like size of mat:
		if (numOfProcesses != numOfPoints) //every process must have a point
		{
			printf("\nInvalid number of processes: %d instead of %d\n", numOfProcesses, numOfPoints);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	MPI_Bcast(&rowSize, 1, MPI_INT, ROOT, MPI_COMM_WORLD); //broadcast the number of points in a row(needed for fuctions)

	MPI_Scatter(mat, 1, PointType, &myPoint, 1, PointType, ROOT, cartComm); //send a point to each proccess.

	if (myId == ROOT)
	{
		printf("\nMatrix as read from file:\n");
		printPointsMatrix(mat, rowSize); //prints matrix before sorting 
		printf("\nAs array:\n");
		printAsArray(mat, rowSize);
	}

	shearSort(coor, &myPoint, PointType, cartComm, rowSize);

	MPI_Gather(&myPoint, 1, PointType, mat, 1, PointType, ROOT, cartComm);

	if (myId == ROOT)
	{
		printf("\nMatrix after sorting:\n");
		printPointsMatrix(mat, rowSize); //prints matrix before sorting 
		printf("\nAs array:\n");
		printAsArray(mat, rowSize);
		free(mat);
	}

	MPI_Finalize();
	return 0;
}

//------------------------------------functions:------------------------------------//
MPI_Comm createCartComm(int rank, int numOfProcesses, int *coords)
{
	MPI_Comm cartComm;
	int nDims = N_DIMS;
	int dim[N_DIMS], period[N_DIMS], reorder;
	int n = (int)sqrt(numOfProcesses);
	dim[0] = n;			// num of columns
	dim[1] = n;			// num of rows
	period[0] = 0;		// cols are not cyclic
	period[1] = 0;		// rows are not cyclic
	reorder = 1;		// allows changing the order of processes ids

	MPI_Cart_create(MPI_COMM_WORLD, nDims, dim, period, reorder, &cartComm);
	MPI_Cart_coords(cartComm, rank, nDims, coords);

	return cartComm;
}

void loadPointsMatrixFromFile(Point *mat, FILE *fp, int numOfPoints)
{
	int i, j;
	int n = (int)sqrt(numOfPoints);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			fscanf_s(fp, "%lf", &mat[i * (int)n + j].x);
			fscanf_s(fp, "%lf", &mat[i * (int)n + j].y);
			fscanf_s(fp, "%lf", &mat[i * (int)n + j].z);
		}
	}
}

double calculateDistanceOfPoint(Point* p1)
{
	double x, y, z;
	x = p1->x;
	y = p1->y;
	z = p1->z;
	return sqrt(x * x + y * y + z * z);
}

void printPointsMatrix(Point* mat, int pointsInRow)
{
	int i, j, n = pointsInRow;

	//prints as a matrix. when sorted it'll be as a snake
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			printf("(%.2f,%.2f,%.2f)\t", mat[i * n + j].x, mat[i * n + j].y, mat[i * n + j].z);
		printf("\n");
	}
}

void printAsArray(Point* mat, int pointsInRow)
{
	int i, j, n = pointsInRow;
	for (i = 0; i < n; i++)
	{
		if (i % 2 == 0)
		{
			for (j = 0; j < n; j++)
			{
				printf("(%.2f,%.2f,%.2f)  ", mat[i * n + j].x, mat[i * n + j].y, mat[i * n + j].z);
			}
		}
		else
		{
			for (j = n - 1; j >= 0; j--)
			{
				printf("(%.2f,%.2f,%.2f)  ", mat[i * n + j].x, mat[i * n + j].y, mat[i * n + j].z);
			}
		}
	}
	printf("\n");
}

void createPointType(MPI_Datatype* PointType)
{
	int blockLengths[POINT_NUM_OF_ATTRIBUTES] = { 1, 1, 1 }; //three attributes in point, each the length of 1
	MPI_Aint disp[POINT_NUM_OF_ATTRIBUTES];
	MPI_Datatype types[POINT_NUM_OF_ATTRIBUTES] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE }; //each attribute is double

	disp[0] = offsetof(Point, x);
	disp[1] = offsetof(Point, y);
	disp[2] = offsetof(Point, z);

	MPI_Type_create_struct(POINT_NUM_OF_ATTRIBUTES, blockLengths, disp, types, PointType);
	MPI_Type_commit(PointType);
}

void shearSort(int* coor, Point* myPoint, MPI_Datatype PointType, MPI_Comm cartComm, int pointsInRow)
{
	const int TAG = 0;
	const int DISP = 1;

	int rank;
	int direction, j;
	int n = pointsInRow * pointsInRow;
	int iterations = (int)ceil(log2(n) + 1); //number of iterations in shear sort, rounded up
	MPI_Status status;
	int source, dest;
	Point hisPoint;

	MPI_Comm_rank(cartComm, &rank);

	for (direction = 0; direction < iterations; direction++) //loop for log2(n)+1
	{
		if (direction % 2 == 0) //sort rows;
		{
			for (j = 0; j < pointsInRow; j++)
			{
				MPI_Cart_shift(cartComm, ROW_DIR, DISP, &source, &dest); //SORTING ROWS
				if (j % 2 == coor[ROW_DIR] % 2)
				{
					if (dest != -1)
					{
						MPI_Send(myPoint, 1, PointType, dest, TAG, cartComm);
						MPI_Recv(myPoint, 1, PointType, dest, TAG, cartComm, &status);
					}
				}
				else
				{
					if (source != -1)
					{
						MPI_Recv(&hisPoint, 1, PointType, source, TAG, cartComm, &status);
						if (coor[COL_DIR] % 2 == 0) //if in even row, sort ascending
						{
							if (calculateDistanceOfPoint(&hisPoint) > calculateDistanceOfPoint(myPoint))
							{
								Point temp = hisPoint;
								hisPoint = *myPoint;
								*myPoint = temp;
							}
						}
						else //if in odd row, sort descending
						{
							if (calculateDistanceOfPoint(&hisPoint) < calculateDistanceOfPoint(myPoint))
							{
								Point temp = hisPoint;
								hisPoint = *myPoint;
								*myPoint = temp;;
							}
						}
						MPI_Send(&hisPoint, 1, PointType, source, TAG, cartComm);
					}
				}
			}
		}
		else //sort columns:
		{
			for (j = 0; j < pointsInRow; j++)
			{
				MPI_Cart_shift(cartComm, COL_DIR, DISP, &source, &dest);

				if (j % 2 == coor[COL_DIR] % 2)
				{
					if (dest != -1)
					{
						MPI_Send(myPoint, 1, PointType, dest, TAG, cartComm);
						MPI_Recv(myPoint, 1, PointType, dest, TAG, cartComm, &status);
					}
				}
				else
				{
					if (source != -1)
					{
						MPI_Recv(&hisPoint, 1, PointType, source, TAG, cartComm, &status);
						if (calculateDistanceOfPoint(&hisPoint) > calculateDistanceOfPoint(myPoint))
						{
							Point temp = hisPoint;
							hisPoint = *myPoint;
							*myPoint = temp;
						}
						MPI_Send(&hisPoint, 1, PointType, source, TAG, cartComm);
					}
				}
			}
		}
		MPI_Barrier(cartComm); //after sorting rows/cols wait for everyone else to finish

	}
}