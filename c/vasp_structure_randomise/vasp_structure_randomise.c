/**
 *	
 *	VASP STRUCTURE RANDOMISE
 * AUTHOR: Artur Tamm <arturt@ut.ee>
 * 
 * This program is free software: you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your option) 
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with 
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include <kdtree.h>

//Because mac os is stupid 
#ifdef MAC
/* ----------------------------- */
/*        T Y P E D E F S        */
/* ----------------------------- */
#define TRUE 1
#define FALSE 0 

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 * getline() -  get line from input stream and return EOF status
 * Derived from the NCSA HTTP util.c code.
 * REturn data in preallocated buffer of max size [0:n-1].
 * Return value of feof(f).
 * Had to redo this function because getline allocates memory
*/
ssize_t getline(char** s, size_t* n, FILE *f)
{/*getline*/
	if(*s == NULL || *n == 0)
	{
		*n = 200;
		*s = malloc((*n) * sizeof(char));
	}
	
	register int nMinus1= (*n-1), i= 0;
	
	while(TRUE)
	{
		register int ch= (char)getc(f);
		
		if(ch == '\r')
		ch= getc(f);
 
		if((ch == '\n') || (ch == EOF) || (i == nMinus1))
		{
			(*s)[i]= ch;
			(*s)[i+1] = '\0';
			return i+1;
		}
		else (*s)[i] = ch;
		
		++i;
		}
}/*getline*/ 
#endif

#define TEST_RESULT(X) 						\
	{												\
		if(X == -1)								\
		{											\
			printf("=====FAILED=====\n");	\
			return -1;							\
		}											\
	}												\

#define TEST_RESULT_MAIN(X,Y)				\
	{												\
		if(X == -1)								\
		{											\
			printf("=====FAILED=====\n");	\
			printf(Y);							\
			printf("=====FAILED=====\n");	\
			cleanup();							\
			return -1;							\
		}											\
	}												\

// atom coordinates
double* atom_x = NULL;
double* atom_y = NULL;
double* atom_z = NULL;

// atom types
int* atom_type = NULL;
int* atom_index = NULL;

// atom neighbour list
int** atom_neighs = NULL;
int*** atom_neigh_index = NULL;
double*** atom_neigh_distance = NULL;

// atom number
int atoms = 0;

// how many types do we have
int atom_types = 0;

// how mane atoms of each type do we have
int* atoms_of_type = NULL;
int ** atom_index_of_type = NULL;

// string consisting of element names (This we don't use only print out in the final result)
// char** atom_names;
char* atom_names = NULL;

// system name
char* system_name = NULL;

// lattice vectors
double vec_a[3] = {0., 0., 0.};
double vec_b[3] = {0., 0., 0.};
double vec_c[3] = {0., 0., 0.};

// steps
int steps = 1000;

// tolerance
double mean_tol = 0.1;

// mirror cells
int mirror_a = 1;
int mirror_b = 1;
int mirror_c = 1;

// up to which neighbour should we look
int neighs = 3;
double* atom_ranges;

// shell radius
double shell = -0.1;

// scaling constant
double scale = 0.;

// kdtree pointer
void* kd = NULL;

// input file name
char* input_file = NULL;

// output file name
char* output_file = NULL;

// find minimal of three numbers
double min3d(double a, double b, double c)
{
	if(a<=b && a<=c) return a;
	else if(b<=a && b<=c) return b;
	else if(c<=a && c<=b) return c;

	return 0.;
}

// read POSCAR file
int read_poscar()
{
	printf("Reading POSCAR file\n=====BEGIN=====\n");

	// read the POSCAR file
	FILE* fd;

	if(input_file == NULL) 
	{
		input_file = malloc(10*sizeof(char));
		strcpy(input_file, "POSCAR");
	}

	fd = fopen(input_file,"r");

	char* line = NULL;
	size_t length;
	ssize_t chars_read;
	int status;

	int i,j,k;

	// read system name
	system_name = NULL;
	length = 0;

	chars_read = getline(&system_name, &length, fd);
	
	TEST_RESULT(chars_read);
	printf("SYSTEM = %s", system_name);

	// read scaling factor
	length = 0;
	chars_read = 0;
	chars_read = getline(&line, &length, fd);
	
	TEST_RESULT(chars_read);
	status = sscanf(line, "%lf", &scale);
	TEST_RESULT(status);
	printf("SCALE  = %10.6f\n", scale);
	
	// read lattice vectors
	chars_read = getline(&line, &length, fd);
	TEST_RESULT(chars_read);
	status = sscanf(line, "%lf %lf %lf", &vec_a[0], &vec_a[1], &vec_a[2]);
	TEST_RESULT(status);
	printf("VECA   = %10.6f %10.6f %10.6f\n", vec_a[0], vec_a[1], vec_a[2]);

	chars_read = getline(&line, &length, fd);
	TEST_RESULT(chars_read);
	status = sscanf(line, "%lf %lf %lf", &vec_b[0], &vec_b[1], &vec_b[2]);
	TEST_RESULT(status);
	printf("VECB   = %10.6f %10.6f %10.6f\n", vec_b[0], vec_b[1], vec_b[2]);

	chars_read = getline(&line, &length, fd);
	TEST_RESULT(chars_read);
	status = sscanf(line, "%lf %lf %lf", &vec_c[0], &vec_c[1], &vec_c[2]);
	TEST_RESULT(status);
	printf("VECC   = %10.6f %10.6f %10.6f\n", vec_c[0], vec_c[1], vec_c[2]);

	// read atom elements
	length = 0;
	chars_read = getline(&atom_names, &length, fd);
	TEST_RESULT(chars_read);
	printf("TYPES  = %s", atom_names);
	
	// read atom numbers of one type
	free(line);
	line = NULL;
	length = 0;
	chars_read = getline(&line, &length, fd);
	
	char* token;
	atom_types = 0;
	
	token = strtok(line, " ");
	if(token == NULL) {
		TEST_RESULT(-1);
	}
	else
	{
		atom_types++;
		atoms_of_type = realloc(atoms_of_type, atom_types*sizeof(int));
		status = sscanf(token, "%d",&atoms_of_type[atom_types-1]);
		if(status == -1) atoms_of_type = realloc(atoms_of_type, (--atom_types)*sizeof(int));

		while(1)
		{
			token = strtok(NULL, " ");
			if(token == NULL) break;

			atom_types++;
			atoms_of_type = realloc(atoms_of_type, atom_types*sizeof(int));
			status = sscanf(token, "%d",&atoms_of_type[atom_types-1]);
			if(status == -1) atoms_of_type = realloc(atoms_of_type, (--atom_types)*sizeof(int));
		}
	}

	printf("ATOMS  =");
	for(i = 0; i < atom_types; i++)
	{
		printf(" %d", atoms_of_type[i]);
		atoms += atoms_of_type[i];
	}

	printf("\nTOTAL  = %d\n", atoms);

	// read coordinate system (not used)
	chars_read = getline(&line, &length, fd);
	TEST_RESULT(chars_read);
	printf("COSYS  = %s", line);

	// initialize arrays
	atom_x = malloc(atoms*sizeof(double));
	atom_y = malloc(atoms*sizeof(double));
	atom_z = malloc(atoms*sizeof(double));
	atom_type = malloc(atoms*sizeof(int));
	atom_index = malloc(atoms*sizeof(int));
	atom_neighs = malloc(atoms*sizeof(int));
	
	atom_index_of_type = malloc(atom_types*sizeof(int*));
	
	for(i = 0; i < atom_types; i++)
	{
		atom_index_of_type[i] = malloc(atoms_of_type[i]*sizeof(int));
	}

	for(i = 0, k = 0; i < atom_types; i++) 
	{
		for(j = 0; j < atoms_of_type[i]; j++, k++)
		{
			chars_read = getline(&line, &length, fd);
			TEST_RESULT(chars_read);
			status = sscanf(line, "%lf %lf %lf", &atom_x[k], &atom_y[k], &atom_z[k]);
			TEST_RESULT(status);

			atom_type[k] = i;
			atom_index[k] = k;

			printf("ATOM   = %4d %2d %10.6f %10.6f %10.6f\n", k, i, atom_x[k], atom_y[k], atom_z[k]);
		}
	}

	for(i = 0; i < atom_types; i++)
	{
		k = 0;

		for(j = 0; j < atoms; j++)
		{
			if(atom_type[j] == i)
			{
				atom_index_of_type[i][k] = atom_index[j];
				k++;
			}
		}

		if(k < atoms_of_type[i]) 
			printf("Warning: Some atoms are missing from list: type %d\n", i);
	}

	if(line != NULL) free(line);
	fclose(fd);
	printf("======END======\n");

	return 0;
}

// write result as poscar file
int write_poscar()
{
	int i,j,k;

	printf("Writing POSCAR file\n");
	printf("=====BEGIN=====\n");

	FILE *fd;

	if(output_file == NULL)
	{
		output_file = malloc(12*sizeof(char));
		strcpy(output_file, "POSCAR_NEW");
	}

	fd = fopen(output_file, "w");

	// write system name
	printf("Writing system name\n");
	fprintf(fd, "%s", system_name);

	// write scaling factor
	printf("Writing scaling factor\n");
	fprintf(fd, "   %-10.6f\n", scale);

	// write lattice vectors
	printf("Writing lattice vectors\n");
	fprintf(fd, "%-10.6f %-10.6f %-10.6f\n", vec_a[0], vec_a[1], vec_a[2]);
	fprintf(fd, "%-10.6f %-10.6f %-10.6f\n", vec_b[0], vec_b[1], vec_b[2]);
	fprintf(fd, "%-10.6f %-10.6f %-10.6f\n", vec_c[0], vec_c[1], vec_c[2]);

	// write atom types
	printf("Writing atom types\n");
	fprintf(fd, "%s", atom_names);

	// write atom type numbers
	printf("Writing atom type numbers\n");
	fprintf(fd, "   ");
	for(i = 0; i < atom_types; i++)
	{
		fprintf(fd, "%-4d ", atoms_of_type[i]);
	}
	fprintf(fd, "\n");

	// write coordinate system
	fprintf(fd, "Direct\n");

	// write atom coordinates
	printf("Writing atom coordinates\n");

	for(i = 0; i < atom_types; i++)
	{
		for(j = 0; j < atoms_of_type[i]; j++)
		{
			fprintf(fd, "%-10.6f   %-10.6f   %-10.6f\n", \
				atom_x[atom_index_of_type[i][j]], \
				atom_y[atom_index_of_type[i][j]], \
				atom_z[atom_index_of_type[i][j]]);
		}
	}

	fclose(fd);


	printf("======END======\n");
	return 0;
}

// find closest neighbours and create neighbour lists
int find_neighs()
{
	int i,j,k,l,m;
	double x,y,z;
	int status;
	
	printf("Finding neighbours\n");
	
	printf("=====BEGIN=====\n");

	//find lattice vector lengths
	double vec_a_length = sqrt(vec_a[0]*vec_a[0] + vec_a[1]*vec_a[1] + vec_a[2]*vec_a[2]);
	double vec_b_length = sqrt(vec_b[0]*vec_b[0] + vec_b[1]*vec_b[1] + vec_b[2]*vec_b[2]);
	double vec_c_length = sqrt(vec_c[0]*vec_c[0] + vec_c[1]*vec_c[1] + vec_c[2]*vec_c[2]);

	double min_length = min3d(vec_a_length, vec_b_length, vec_c_length);

	//create kdtree
	void *kd;
	kd = kd_create(3);

	//FILE *fd;
	//fd = fopen("TEMP.xyz","w");
	
	m=0;
	//insert atoms
	for(i = 0; i < atoms; i++)
	{
		for(j = -mirror_a; j <= mirror_a; j++)
		{
			for(k = -mirror_b; k <= mirror_b; k++)
			{
				for(l = -mirror_c; l <= mirror_c; l++)
				{
					x = scale*(atom_x[i]*vec_a[0] + j*vec_a[0] + \
						atom_y[i]*vec_b[0] + k*vec_b[0] + \
						atom_z[i]*vec_c[0] + l*vec_c[0]);

					y = scale*(atom_x[i]*vec_a[1] + j*vec_a[1] + \
						atom_y[i]*vec_b[1] + k*vec_b[1] + \
						atom_z[i]*vec_c[1] + l*vec_c[1]);

					z = scale*(atom_x[i]*vec_a[2] + j*vec_a[2] + \
						atom_y[i]*vec_b[2] + k*vec_b[2] + \
						atom_z[i]*vec_c[2] + l*vec_c[2]);

					kd_insert3(kd, x, y, z, &atom_index[i]);

					m++;
	//				fprintf(fd, "%d %10.6f %10.6f %10.6f\n", \
	//					atom_type[atom_index[i]], x, y, z);
				}
			}
		}
	}

	//printf("We have %d atoms in MiRROR WORLD\n",m);

	//fclose(fd);

	void *kdresult;
	double range = 0.;
	double range_step = 0.101; //in Angstroms
	double length = 0.0;

	// allocate the arrays
	atom_neighs = malloc(neighs*sizeof(int*));
	atom_neigh_index = malloc(neighs*sizeof(int**));
	atom_neigh_distance = malloc(neighs*sizeof(double**));
	atom_ranges = malloc(neighs*sizeof(double*));

	for(i = 0; i < neighs; i++)
	{
		atom_neighs[i] = malloc(atoms*sizeof(int));
		atom_neigh_index[i] = malloc(atoms*sizeof(int*));
		atom_neigh_distance[i] = malloc(atoms*sizeof(double*));
	}
	
	i=0;
	j=0;
	k=0;

	x = scale*(atom_x[0]*vec_a[0] + atom_y[0]*vec_b[0] + atom_z[0]*vec_c[0]);
	y = scale*(atom_x[0]*vec_a[1] + atom_y[0]*vec_b[1] + atom_z[0]*vec_c[1]);
	z = scale*(atom_x[0]*vec_a[2] + atom_y[0]*vec_b[2] + atom_z[0]*vec_c[2]);
		
	// find the requested neighbour distances
	// this one is not nice when we have distance calculation errors similar to 
	// range_step
	// RETHINK TIHS ONE (give neighbour list numbers)
	// IN THIS CASE WE COULD ALSO HAVE SHELL RADIUSES IN THE SAME OPTION
	if(shell > 0.0) 
	{
		printf("NUMBER OF SHELLS IS %2d WITH STEPS OF SIZE %6.2f\n", neighs, shell);

		for(i=1; i<=neighs; i++)
		{
			atom_ranges[i-1] = shell*i;
		}
	}
	else
	{
		while(1)
		{
			kdresult = kd_nearest_range3(kd, x, y, z, range);
			j = kd_res_size(kdresult);
			kd_res_free(kdresult);
	
			if(i != j) {
				printf("NUMBER OF %2d NEIGHBOURS = %4d UPTO DISTANCE %10.6f \n",k,j,range);
				i = j;
				if(k != 0) atom_ranges[k-1] = range;
	
				k++;
			}
			else
			{
				range += range_step;
			}
	
			if(k > neighs) break;
		}
	}

	// now create the neighbour lists for all the atoms
	// THIS MEGA THING WORKS
	for(i = 0; i < neighs; i++)
	{
		for(j = 0; j < atoms; j++)
		{
			x = scale*(atom_x[j]*vec_a[0] + atom_y[j]*vec_b[0] + atom_z[j]*vec_c[0]);
			y = scale*(atom_x[j]*vec_a[1] + atom_y[j]*vec_b[1] + atom_z[j]*vec_c[1]);
			z = scale*(atom_x[j]*vec_a[2] + atom_y[j]*vec_b[2] + atom_z[j]*vec_c[2]);

			kdresult = kd_nearest_range3(kd, x, y, z, atom_ranges[i]);
			
			atom_neighs[i][j] = kd_res_size(kdresult) - 1;
			
			//printf("ATOM %d LEVEL %d NEIGHS %d POINT %10.6f %10.6f %10.6f\n" \
			//	,j ,i, atom_neighs[i][j], x, y, z);

			if(i != 0)
			{
				for(k = i - 1; k >= 0; k--)
					atom_neighs[i][j] -= atom_neighs[k][j];
			}

			atom_neigh_index[i][j] = malloc(atom_neighs[i][j]*sizeof(int));
			atom_neigh_distance[i][j] = malloc(atom_neighs[i][j]*sizeof(double));

			int* index;
			double xx, yy, zz;
			l = 0;

			while(1)
			{
				xx = 1.; yy = 1.; zz = 1.;
				kd_res_item3(kdresult, &xx, &yy, &zz);
				index = kd_res_item_data(kdresult);
				length = sqrt((xx - x)*(xx - x) + (yy - y)*(yy - y) + (zz - z)*(zz - z));

				if((i != 0 && length > atom_ranges[i - 1]) || (i == 0 && length > 0.0))
				{
					atom_neigh_index[i][j][l] = *index;
					atom_neigh_distance[i][j][l] = length;
					l++;
				}

				if(kd_res_next(kdresult) == 0) break;
			}

			kd_res_free(kdresult);
		}
	}

	printf("======END======\n");

	return 0;
}

// find the correlation of neigh'th neighbour
void find_corr(int neigh, double* corr)
{
	double atom_corr = 0.;
	int i,j,k;

	for(i = 0; i < atom_types; i++)
	{
		corr[i] = 0.;
	}

	for(i = 0; i < atoms; i++)
	{
		atom_corr = 0.;

		for(j = 0; j < atom_neighs[neigh][i]; j++)
		{
			if(atom_type[atom_neigh_index[neigh][i][j]] == atom_type[i])
			{
				atom_corr += 1.;
			}
			else 
			{
				atom_corr -= 1.;
			}
		}

		corr[atom_type[i]] += atom_corr / (double) atom_neighs[neigh][i];
		//printf("%d %d %d\n",neigh, i, atom_neighs[neigh][i]);
	}
	
	for(i = 0; i < atom_types; i++)
	{
		corr[i] = corr[i] / (double) atoms_of_type[i];
	}
}

int iterate()
{
	int best_type[atoms];
	int** best_index_of_type;
	double best_corr[neighs][atom_types];
	double best_mean;
	
	double current_corr[neighs][atom_types];
	double current_mean;

	double theor_corr[atom_types];
	int i,j,k,l,m,n,p;
	int status;

	printf("Starting iterations\n");
	printf("=====BEGIN=====\n");

	// here we could improve the program by finding different correlations 
	// for differet types

	best_index_of_type = malloc(atom_types*sizeof(int*));

	for(i = 0; i < atom_types; i++)
	{
		best_index_of_type[i] = malloc(atoms_of_type[i]*sizeof(int));
	}

	// find theoretical correlation
	printf("THEORETICAL CORRELATION = ");

	for(i = 0; i < atom_types; i++)
	{
		theor_corr[i] = (2.0 * (double) atoms_of_type[i] / ((double) atoms) - 1.0);
		//theor_corr[i] *= theor_corr[i];
		printf("%10.6f ", theor_corr[i]);
	}

	printf("\n");

	// find correlation of the current system
	printf("CURRENT CORRELATION = \n");
	
	best_mean = 0.;

	for(i = 0; i < neighs; i++)
	{
		find_corr(i, best_corr[i]);
		printf("NEIGH = %d ", i);

		for(j = 0; j < atom_types; j++)
		{
			best_mean += (theor_corr[j] - best_corr[i][j]) * \
				(theor_corr[j] - best_corr[i][j]);
			printf("%10.6f ", best_corr[i][j]);
		}

		printf("\n");
	}

	// best_mean = best_mean / atom_types / neighs; // this is not necessary

	// print mean square
	printf("CURRENT MEAN = %10.6f\n", best_mean);

	// copy current atom types to the best one
	memcpy(best_type, atom_type, atoms*sizeof(int));
	
	// copy the vectors
	for(i = 0; i < atom_types; i++)
	{
		memcpy(best_index_of_type[i], atom_index_of_type[i], atoms_of_type[i]*sizeof(int));
	}

	srandom((unsigned int) time(NULL));

	// this could be parallelised
	while(i < steps)
	{
		l = random()%atom_types;
		j = random()%atoms_of_type[l];
		
		m = l;
		while(m == l)
			m = random()%atom_types;
		
		k = random()%atoms_of_type[m];
		
		//printf("%d-%d\n", l, m);

		n = atom_type[atom_index_of_type[l][j]];
		atom_type[atom_index_of_type[l][j]] = atom_type[atom_index_of_type[m][k]];
		atom_type[atom_index_of_type[m][k]] = n;

		n = atom_index_of_type[l][j];
		atom_index_of_type[l][j] = atom_index_of_type[m][k];
		atom_index_of_type[m][k] = n;

		for(j = 0; j < neighs; j++)
		{
			find_corr(j, current_corr[j]);
		}

		current_mean = 0.;

		for(j = 0; j < neighs; j++)
		{
			for(k = 0; k < atom_types; k++)
			{
				current_mean += (theor_corr[k] - current_corr[j][k]) * \
					(theor_corr[k] - current_corr[j][k]);
			}
		}

		if(current_mean < best_mean)
		{
			memcpy(best_type, atom_type, atoms*sizeof(int));
			
			for(k = 0; k < atom_types; k++)
			{
				memcpy(best_index_of_type[k], atom_index_of_type[k], atoms_of_type[k]*sizeof(int));
			}

			best_mean = current_mean;
		}

		i++;

		if(best_mean < mean_tol) break;
	}

	if(i == steps) status = 1;
	else status = 0;

	// copy best result to output
	memcpy(atom_type, best_type, atoms*sizeof(int));

	for(i = 0; i < atom_types; i++)
	{
		memcpy(atom_index_of_type[i], best_index_of_type[i], atoms_of_type[i]*sizeof(int));
	}

	// find correlation of the current system
	printf("FINAL CORRELATION = \n");

	best_mean = 0.;

	for(i = 0; i < neighs; i++)
	{
		find_corr(i, best_corr[i]);
		printf("NEIGH = %d ", i);

		for(j = 0; j < atom_types; j++)
		{
			best_mean += (theor_corr[j] - best_corr[i][j]) * \
				(theor_corr[j] - best_corr[i][j]);
			printf("%10.6f ", best_corr[i][j]);
		}

		printf("\n");
	}

	// print mean square
	printf("CURRENT MEAN = %10.6f\n", best_mean);
	
	printf("======END======\n");

	// release allocated memory
	for(i = 0; i < atom_types; i++)
	{
		if(best_index_of_type[i] != NULL) free(best_index_of_type[i]);
	}

	free(best_index_of_type);

	return status;
}

// free used memory
void cleanup()
{
	int i,j;
	
	if(atom_x != NULL) free(atom_x);
	if(atom_y != NULL) free(atom_y);
	if(atom_z != NULL) free(atom_z);
	if(atom_type != NULL) free(atom_type);
	if(atom_index != NULL) free(atom_index);
	if(atoms_of_type != NULL) free(atoms_of_type);

	if(output_file != NULL) free(output_file);
	if(input_file != NULL) free(input_file);

	if(atom_neighs != NULL)
	{
		for(i = 0; i < neighs; i++) 
		{
			if(atom_neighs[i] != NULL) free(atom_neighs[i]);
		}

		free(atom_neighs);
	}

	if(atom_neigh_index != NULL) 
	{
		for(i = 0; i < neighs; i++) 
		{
			if(atom_neigh_index[i] != NULL)
			{
				for(j = 0; j < atoms; j++) 
					if(atom_neigh_index[i][j] != NULL) free(atom_neigh_index[i][j]);

				free(atom_neigh_index[i]);
			}
		}

		free(atom_neigh_index);
	}

	if(atom_neigh_distance != NULL) 
	{
		for(i = 0; i < neighs; i++) 
		{
			if(atom_neigh_distance[i] != NULL)
			{
				for(j = 0; j < atoms; j++)
					if(atom_neigh_distance[i][j] != NULL) free(atom_neigh_distance[i][j]);

				free(atom_neigh_distance[i]);
			}
		}

		free(atom_neigh_distance);
	}

	if(atom_index_of_type != NULL)
	{
		for(i = 0; i < atom_types; i++)
			if(atom_index_of_type[i] != NULL) free(atom_index_of_type[i]);

		free(atom_index_of_type);
	}

	if(atom_names != NULL) free(atom_names);
	if(system_name != NULL) free(system_name);

	if(kd != NULL) kd_free(kd);
}

int main(int argc, char** argv)
{
	// this is just for status of some functions
	int status;
	int opt;

	status = 0;

	// put option parsing here
	while((opt = getopt(argc, argv, "a:b:c:hs:2:t:i:o:g:")) != -1)
	{
		switch(opt) 
		{
			case 'a':
				mirror_a = atoi(optarg);
				if(mirror_a < 1) mirror_a = 1;
				break;
			case 'b':
				mirror_b = atoi(optarg);
				if(mirror_b < 1) mirror_b = 1;
				break;
			case 'c':
				mirror_c = atoi(optarg);
				if(mirror_c < 1) mirror_c = 1;
				break;
			case '2':
				neighs = atoi(optarg);
				if(neighs < 1) neighs = 1;
				break;
			case 'g':
				shell = atof(optarg);
				if(shell < 0.0) shell = 0.1;
				break;
			case 's':
				steps = atoi(optarg);
				if(steps < 1) steps = 1;
				break;
			case 't':
				mean_tol = atof(optarg);
				if(mean_tol < 0.) mean_tol = 0.;
				break;
			case 'i':
				input_file = malloc(sizeof(char)*(strlen(optarg) + 1));
				strcpy(input_file, optarg);
				break;
			case 'o':
				output_file = malloc(sizeof(char)*(strlen(optarg) + 1));
				strcpy(output_file, optarg);
				break;
			case 'h':
				printf("Vasp structure randomise program by Artur Tamm <arturt@ut.ee>\n");
				printf("Program help:\n");
				printf("\t-h\t\t\tdisplay this help and exit\n");
				printf("\t-2 [int]\t\tcalculate correlation for upto this neighbour/shell: default 3\n");
				printf("\t-g [float]\t\tcalculate correlation in a shell with this radius: default use nearest neighbours\n");
				printf("\t-a [int]\t\tmake [-a .. a] copies in lattice vector a direction; default 1\n");
				printf("\t-b [int]\t\tmake [-b .. b] copies in lattice vector b direction; default 1\n");
				printf("\t-c [int]\t\tmake [-c .. c] copies in lattice vector c direction; default 1\n");
				printf("\t-s [int]\t\tdo at most this many steps; default 1000\n");
				printf("\t-t [float]\t\tend iteration when this tolerance is met: default 0.1\n");
				printf("\t-i [string]\t\tvasp structure input file name; default POSCAR\n");
				printf("\t-o [string]\t\tvasp structure output file name; default POSCAR_NEW\n");

				printf("\n");
				status = -1;
				break;
			default:
				status = -1;
		}

		if(status == -1) break;
	}

	// print usage and exit
	if(status == -1) 
	{
		printf("Usage:\n vasp_structure_randomise -2 [int] -g [float] -a [int] -b [int] -c [int] -s [int] -t [float] -i [string] -o [string]\n");
		return 0;
	}

	// start reading the POSCAR file 
	status = read_poscar();
	TEST_RESULT_MAIN(status,"read_poscar()\n");

	// create neighbour lists
	status = find_neighs();
	TEST_RESULT_MAIN(status, "find_neighs()\n");
   
	// start iterating
	status = iterate();
	if(status == 0)
	{
		printf("Reached required tolerance\n");
	}
	else if(status == 1)
	{
		printf("Reached maximum number of iterations\n");
	}

	// write result to output
	status = write_poscar();

	// clean everything up
	cleanup();
}

