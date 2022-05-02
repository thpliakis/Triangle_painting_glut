#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

// #define MAX(x, y) (((x) > (y)) ? (x) : (y))
// #define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define ALPHA 0.5
typedef struct 
{
	GLint v0[2];
	GLint v1[2];
	GLint v2[2];
}triangle;

typedef struct 
{
	GLfloat c0[3];
	GLfloat c1[3];
	GLfloat c2[3];
}color;

typedef struct
{
	int verts2d[4999][2];
	float vcolors[4999][3];
	int faces[10000][3];
	float depth[4999][1];
}data;

data d;

typedef struct 
{
	GLint x;
	GLint y;
	GLfloat rgb[3];
}point;


// float interpolate_color(point p1, point p2, point p){
// 	float r1  = sqrt(pow(p2.x - p.x,2) + pow(p2.y - p.y,2));
// 	float r2  = sqrt(pow(p.x - p1.x,2) + pow(p.y - p1.y,2));
// 	float percent = r1 / (r1+r2);
// 	for(int i=0; i<3; i++){
// 		p.rgb[i] = percent * p1.rgb[i] + (1 - percent) * p2.rgb[i];
// 	}
// }

void interpolate_color(GLint *p1, GLint *p2, GLint *p, GLfloat *c1,GLfloat *c2, GLfloat *c){
	GLfloat r1  = sqrt(pow(p2[0] - p[0],2) + pow(p2[1] - p[1],2));
	GLfloat r2  = sqrt(pow(p[0] - p1[0],2) + pow(p[1] - p1[1],2));
	GLfloat percent = r1 / (r1+r2);
	for(int i=0; i<3; i++){
		c[i] = percent * c1[i] + (1 - percent) * c2[i];
	}
}

void read_data_float( int M, int N, float m[M][N],  FILE *file){

	char line[250];

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			if(fgets(line, sizeof(line), file)){
				m[i][j] = atof(line);
			}else {
				break;
			}
		}
	}
}

void read_data_int( int M, int N, int m[M][N],  FILE *file){

	char line[250];

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			if(fgets(line, sizeof(line), file)){
				m[i][j] = (int)atof(line);
			}else {
				break;
			}
		}
	}
}

void swap_int(int *a, int k, int l) {
	int temp = a[k];
	a[k] = a[l];
	a[l] = temp;
}

void swap_float(float *a, int k, int l) {
	float temp = a[k];
	a[k] = a[l];
	a[l] = temp;
}


/* partition -- in-place update of elements */
int partition2(float *a,int *o, int n) {
	float pivot = a[n-1];
	int i = 0;
	
	for (int j = 0; j < n - 1; j++){
		if (a[j] >= pivot){ 
		swap_int(o,i,j);
		swap_float(a,i++,j);
		}
	}

	swap_float(a, i, n - 1);
	swap_int(o, i, n - 1);
	return (i);
}

/* qsortseq -- Entry point for QuickSort */
void q_sort(float *a,int *o, int n) {
		
	if (n > 1) {
		int p = partition2(a,o, n);
		q_sort(a,o,p);
		q_sort(&a[p+1],&o[p+1],n-p-1);
	}
}

GLint min(int l, GLint *A){
	GLfloat min = A[0];
	for(int i=1; i<l; i++){
		if(A[i]<min)
			min = A[i];
	}
	return min;
}

GLint max(int l, GLint *A){
	GLfloat max = A[0];
	for(int i=1; i<l; i++){
		if(A[i]>max)
			max = A[i];
	}
	return max;
}

GLint Min_f(GLint a, GLint b){
	if(a>b)
		return b;
	return a;
}
GLint Max_f(GLint a, GLint b){
	if(a<b)
		return b;
	return a;
}

void shade_triangle(triangle *v, color *c, char *shade_t){
	
	GLfloat a =1;// (GLfloat)rand() / (GLfloat)RAND_MAX;

	GLint Ykmin[3], Xkmin[3],Ykmax[3],Xkmax[3],activeSides[3],Min,Max;
	GLint Xmin, Xmax, Ymin , Ymax;

	GLint x1[2],x2[2],x3[2],x4[2], V[3][2];
	GLfloat C1[3],C2[3],C3[3],C4[3],Xk[3],C[3][3];

	GLint p[2];
	GLfloat cp[3];

	// if(!strcmp(shade_t, "flat")){
	// 	printf("%s",shade_t);
	// }
	// Extra variables for shate_t = "gouraud"
	// Here the x coordinate of the Active side is stored
    // Indexes of the active sides, to know which of the 3 are active
	GLint Xact[3],actIndex[2],Xact2[2],indexAct[2], XminLocal,XmaxLocal;

	// Here the average color of a trinagle is stored
	GLfloat cflat[3],m[3], er[3];

	// The interpolated colors are saved here
	GLfloat colors[3], colosXact[2][3], colorsInside[3];

	// In indexmin and indexman the vertices that every side is between are stored  
    GLint indexmin[3] ,indexmax[3];

	GLint vertex1, vertex2, pairs[3][2]; // Temp variables
	int ptemp=0,cross_c=0;

//     GLfloat C[3][3] = {{0.818012,  0.251958,  0.817453}, 
// 					   {0.923060,  0.989394,  0.942833}, 
//  					   {0.401391,  0.168952,  0.269701}};
// ;
//     GLint V[3][2] = {{185, 172},
//                       {355, 350},
//                       {450, 100}};
//     srand(time(0));

    // for(int i = 0; i<3; i++){
    //     for(int j =0; j<2; j++){
    //         printf(" %d ", V[i][j]);
    //     }
    //     printf("\n");
    // }

	for(int i=0; i<3; i++){
		activeSides[i] = 0;
		er[i] = 0;
        Xact[i] = 0;
		Xk[i] = 0;
		cflat[i];
	}

	V[0][0] = v->v0[0];
	V[0][1] = v->v0[1];
	V[1][0] = v->v1[0];
	V[1][1] = v->v1[1];
	V[2][0] = v->v2[0];
	V[2][1] = v->v2[1];
	C[0][0] = c->c0[0];
	C[0][1] = c->c0[1];
	C[0][2] = c->c0[2];
	C[1][0] = c->c1[0];
	C[1][1] = c->c1[1];
	C[1][2] = c->c1[2];
	C[2][0] = c->c2[0];
	C[2][1] = c->c2[1];
	C[2][2] = c->c2[2];
	// For each side we find the min and max X,Y coordinates and calculate the slope m
	for(int i=0; i<3; i++){
		vertex1 = i;
		// printf("vertex1 = %d\n",vertex1);
		vertex2 = (i==2) ? 0 : i+1 ;
		// printf("vertex2 = %d\n",vertex2);
    	pairs[i][0] = vertex1;
    	pairs[i][1] = vertex2;

		m[i]=(GLfloat)(V[vertex2][1]-V[vertex1][1])/ (GLfloat)(V[vertex2][0]-V[vertex1][0]);
		// printf("m[%d] = %f\n",i,m[i]);

		// Xkmin[i] = Min_f(V[vertex1][0],V[vertex2][0]);
		// Xkmax[i] = Max_f(V[vertex1][0],V[vertex2][0]);
		Ykmin[i] = Min_f(V[vertex1][1],V[vertex2][1]);
		if(Ykmin[i] == V[vertex1][1]){
			Xkmin[i] = V[vertex1][0];
		}else{
			Xkmin[i] = V[vertex2][0];
		}
		Ykmax[i] = Max_f(V[vertex1][1],V[vertex2][1]);
		if(Ykmax[i] == V[vertex1][1]){
			Xkmax[i] = V[vertex1][0];
		}else{
			Xkmax[i] = V[vertex2][0];
		}

	}

	Ymin = V[0][1];
	Xmin = V[0][0];
	Ymax = V[0][1];
	Xmax = V[0][0];

	// for(int i=0; i<3; i++){
		// printf("activeSides[%d] = %d\n",i,activeSides[i]);
		// printf("m[%d] = %d\n",i,m[i]);
	// }

	for(int i=0; i<3; i++){
		// printf("Ykmin[%d] = %d\n",i, Ykmin[i]);
		// printf("Xkmin[%d] = %d\n",i, Xkmin[i]);
		if(V[i][1]<Ymin)
			Ymin = V[i][1];
		if(V[i][0]<Xmin)
			Xmin = V[i][0];
		if(V[i][0]>Xmax)
            Xmax = V[i][0];
		if(V[i][1]>Ymax)
			Ymax = V[i][1];
	}
	printf("Ymin = %d\n", Ymin);
	printf("Ymax = %d\n", Ymax);
	printf("Xmin = %d\n", Xmin);
	printf("Xmax = %d\n\n", Xmax);

	for(int i=0; i<3; i++){
		cflat[i] = 0;
		for(int j=0; j<3; j++)
			cflat[i] += C[j][i];
		cflat[i] = cflat[i] / 3;
	}

	// Ymin = min(3,Ykmin);
	// Ymax = max(3,Ykmax);
	// Xmin = min(3,Xkmin);
	// Xmax = max(3,Xkmax);

	if(Ymin == Ymax){
		if((Xkmax[0]==Xmax) && (Xkmax[1]==Xmax) && (Xkmax[2]==Xmax)){ //if triangle is a point
			glBegin(GL_POINTS);
  				glColor3f(cflat[0],cflat[1],cflat[2]);
  		  		glVertex2i(Ymin,Xmax);
			glEnd();
		}else{ //if triangle is a flat line
			for(int x=Xmin; x<=Xmax; x++){
				ptemp=0;
				for(GLint k=0; k<3; k++){
					Min = Min_f(Xkmin[k],Xkmax[k]);
					Max = Max_f(Xkmin[k],Xkmax[k]);
					if(x>= Min && x<Max){
						Xact2[ptemp] = k;
						activeSides[k] = 1;
						ptemp += 1;
					}else if(x != Xmax)
						activeSides[k] = 0;
				}
				if(!strcmp(shade_t, "gouraud")){
					// interpolate_color(t[Xact2[0]],t[Xact2[1]],p);
					p[0] = x;
					p[1] = Ymin;
					interpolate_color(V[Xact2[0]],V[Xact2[1]],p,C[Xact2[0]],C[Xact2[1]],cp);
					glBegin(GL_POINTS);
						glColor3f(cp[0],cp[1],cp[2]);
						glVertex2i(Ymin,x);
					glEnd();
				}else if (!strcmp(shade_t, "flat")){
					glBegin(GL_POINTS);
						glColor3f(cflat[0],cflat[1],cflat[2]);
						glVertex2i(Ymin,x);
					glEnd();
				}
			
			}
		}

	}
	else{
		// printf("%s\n",shade_t);
		for(int k=0; k<3; k++){
			if(Ykmin[k] == Ymin && m[k] !=0){
				activeSides[k] = 1;
				Xk[k] = ((GLfloat)(Ymin - V[k][1]) / (GLfloat) m[k]) + (GLfloat) V[k][0];
			}

			if(Ykmax[k] == Ymin)
				activeSides[k] = 0;
		}
		

		for(int y=Ymin; y<=Ymax; y++){
			// for(int i=0; i<3; i++)
				// printf("activeSides[%d] = %d\n",i,activeSides[i]);

			for(int i=0; i<3; i++)
				Xact[i] = round(Xk[i]);

			for(int i=0; i<3; i++){
                er[i] = er[i] + (Xk[i] - Xact[i]);
                if(er[i] > 1){    
                    Xact[i] += 1;
                    er[i] = 0;
				}
                if(er[i] < -1){
                    Xact[i] -= 1;
                    er[i] = 0;
				}
			}

			if(activeSides[0]==0 && activeSides[1]==0 && activeSides[2]==0)
                continue;;
            if(activeSides[0]==1 && activeSides[1]==1 && activeSides[2]==1)
                continue;
			

            ptemp=0;
			for(int i=0; i<3; i++){
                if(activeSides[i] == 1){
                    indexAct[ptemp] = i;
                    Xact2[ptemp] = Xact[i];
                    ptemp += 1;
				}
			}
			// printf("Xact2[0] = %d\n",Xact2[0]);
			// printf("Xact2[1] = %d\n",Xact2[1]);
			// for(int i=0; i<3; i++)
			// 	printf("activeSides[%d] = %d\n",i,activeSides[i]);

			XmaxLocal = 0;
			XminLocal = 600;
			for(int i=0; i<2; i++){
					if(Xact2[i]>XmaxLocal){
						XmaxLocal = Xact2[i];
					}
					if(Xact2[i]<XminLocal){
						XminLocal = Xact2[i];
					}
			}
			// printf("XmaxLocal = %d\n",XmaxLocal);
			// printf("XminLocal = %d\n",XminLocal);
			

            // x1,x2,x3,x4 are the vertices that our active points lie between and their respective colors C1,C2,C3,C4
			for(int i=0; i<2; i++){
				x1[i] = V[pairs[indexAct[0]][0]][i];
				x2[i] = V[pairs[indexAct[0]][1]][i];
				x3[i] = V[pairs[indexAct[1]][0]][i];
				x4[i] = V[pairs[indexAct[1]][1]][i];
			}
			for(int i=0; i<3; i++){
				C1[i] = C[pairs[indexAct[0]][0]][i];
				C2[i] = C[pairs[indexAct[0]][1]][i];
				C3[i] = C[pairs[indexAct[1]][0]][i];
				C4[i] = C[pairs[indexAct[1]][1]][i];
			}
            // Calculate the colors of the active points
			p[0] = Xact2[0];
			p[1] = Ymin;
            interpolate_color(x1,x2,p,C1,C2,colosXact[0]);
			p[0] = Xact2[1];
            interpolate_color(x3,x4,p,C3,C4,colosXact[1]);

            // temp variable to correnctly keep track of Xact (active points)
            cross_c=0;      
			// printf("XmaxLocal = %d\n",XmaxLocal);
			// printf("XminLocal = %d\n",XminLocal);
            for(GLint x=XminLocal; x<=XmaxLocal; x++){
                
                if(Xact2[0] == Xact2[1]){
                    if(!strcmp(shade_t, "gouraud")){
						glBegin(GL_POINTS);
							glColor3f(colosXact[0][0],colosXact[0][1],colosXact[0][2]);
							glVertex2i(y,x);
						glEnd();
					}
                    else if(!strcmp(shade_t, "flat")){
						glBegin(GL_POINTS);
							glColor3f(cflat[0],cflat[1],cflat[2]);
							glVertex2i(y,x);
						glEnd();
						//break
					}
				}

                if(Xact2[0] == x){
                    if(!strcmp(shade_t, "gouraud")){
						glBegin(GL_POINTS);
							glColor3f(colosXact[0][0],colosXact[0][1],colosXact[0][2]);
							glVertex2i(y,x);
						glEnd();
					}
                    else if(!strcmp(shade_t, "flat")){
						glBegin(GL_POINTS);
							glColor3f(cflat[0],cflat[1],cflat[2]);
							glVertex2i(y,x);
						glEnd();
					}
                    cross_c += 1; // +1 if we meet a Active side
				}
				else if(Xact2[1] == x){
                    if(!strcmp(shade_t, "gouraud")){
						glBegin(GL_POINTS);
							glColor3f(colosXact[1][0],colosXact[1][1],colosXact[1][2]);
							glVertex2i(y,x);
						glEnd();
					}
                    else if(!strcmp(shade_t, "flat")){
						glBegin(GL_POINTS);
							glColor3f(cflat[0],cflat[1],cflat[2]);
							glVertex2i(y,x);
						glEnd();
					}
                    cross_c += 1; // +1 if we meet a Active side
				}
            
                if((cross_c % 2) == 1){ // When cross_c is an odd number we are inside a triangle
                    if(!strcmp(shade_t, "gouraud")){
						p[0] = x; 
						p[1] = y;

						x1[0] = Xact2[0];
						x2[0] = Xact2[1];
						x1[1] = y;
						x2[1] = y;

						for(int i=0; i<3; i++){
							C1[i] = colosXact[0][i];
							C2[i] = colosXact[1][i];
						}

						interpolate_color(x1,x2,p,C1,C2,colorsInside);
						glBegin(GL_POINTS);
							glColor3f(colorsInside[0],colorsInside[1],colorsInside[2]);
							glVertex2i(y,x);
						glEnd();
					}
                    else if(!strcmp(shade_t, "flat")){
						// printf("flat\n");
						glBegin(GL_POINTS);
							glColor3f(cflat[0],cflat[1],cflat[2]);
							glVertex2i(y,x);
						glEnd();
					}
				}
			}

			// Update the active points with the correct gradient
			for(int k=0; k<3; k++){
				if(activeSides[k] == 1){  
					if(isinf(m[k]))
						Xk[k]=Xk[k];
					if(m[k] == 0)
						Xk[k] += 1;
					else
						Xk[k] = Xk[k] + 1/m[k];
				}
			}
			for(int k=0; k<3; k++){
				if(Ykmin[k] == y+1 && m[k] != 0){
					activeSides[k] = 1;          // When the scanline cuts the side of the triangle
					Xk[k] = (y+1 - (GLfloat) V[k][1]) / m[k] + (GLfloat)V[k][0];
				}
				if((Ykmax[k] == y+1) && (y+1 != Ymax))
					activeSides[k] = 0;
			}
		}
	}
	
	// glBegin(GL_TRIANGLES);
	// 	glColor4f(c->c0[0],c->c0[1],c->c0[2],a);
    //     glVertex2i(v->v1[0],v->v1[1]);
	// 	glColor4f(c->c1[0],c->c1[1],c->c1[2],a);
    //     glVertex2i(v->v1[0],v->v2[1]);
	// 	glColor4f(c->c2[0],c->c2[1],c->c2[2],a);
    //     glVertex2i(v->v2[0],v->v2[1]);
    // glEnd();
	printf("painted\n");
}

void render(void)
{
	triangle v;
	color c;
	char shade_t[8] = "flat";

	glEnable( GL_BLEND );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	
	glClearColor(0.5f, 0.5f, 0.5f, 0.5f);	// Grey ackround color
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);      
	int triang_num = sizeof(d.faces) / sizeof(d.faces[0]);

	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0,600.0,600.0,0.0);
	
	for(int t=0; t<triang_num; t++){
	// for(int t=0; t<1; t++){
		for(int i = 0; i<2; i++){
			v.v0[i] = (GLint)d.verts2d[d.faces[t][0]][i];
			v.v1[i] = (GLint)d.verts2d[d.faces[t][1]][i];
			v.v2[i] = (GLint)d.verts2d[d.faces[t][2]][i];
		}
		
		for(int j=0; j<3; j++){	
			c.c0[j] = (GLfloat)d.vcolors[d.faces[t][0]][j];
			c.c1[j] = (GLfloat)d.vcolors[d.faces[t][0]][j];
			c.c2[j] = (GLfloat)d.vcolors[d.faces[t][0]][j];
			
		}		
		shade_triangle(&v, &c, shade_t);
		
	}

glFlush();
glutSwapBuffers();
//glutPostRedisplay();
sleep(1);
    
}

int main(int argc, char** argv)
{	
	char const* const fileName = "verts2d.txt"; /* should check that argc > 1 */
	char const* const fileName2 = "vcolors.txt";
	char const* const fileName3 = "faces.txt";
	char const* const fileName4 = "depth.txt";

	int temp_f[10000][3];

    FILE *file , *file2, *file3, *file4;
	// = fopen(fileName, "r"); /* should check the result */
	
    if ((file = fopen(fileName, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }
	if ((file2 = fopen(fileName2, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }
	if ((file3 = fopen(fileName3, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }
	if ((file4 = fopen(fileName4, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }

	// Read data from txt files.
	read_data_int(4999,2, d.verts2d,file);
	read_data_float(4999, 3, d.vcolors,file2);
	read_data_int(10000, 3, d.faces, file3);
	read_data_float(4999, 1, d.depth, file4);
	
	int triang_num = sizeof(d.faces) / sizeof(d.faces[0]);
	float *d_Median = (float *) malloc(triang_num * sizeof(float));
	for (int l = 0; l < triang_num; l++) d_Median[l] = 0.0;
	int *d_Order = (int *) malloc(triang_num * sizeof(int)); 
	for (int l = 0; l < triang_num; l++) d_Order[l] = l;

	//Calculate the depth of each triangle (the median of the depths of each vertex)
	for(int i=0; i<triang_num;i++){
		for(int j=0; j<3; j++){
			d_Median[i] = d.depth[d.faces[i][j]][0] + d_Median[i]; 
		}
		d_Median[i] = d_Median[i]/3;
	}
	q_sort(d_Median,d_Order,triang_num);

	for(int i=0; i<triang_num;i++)
		for(int j=0; j<3; j++)	
			temp_f[i][j] = d.faces[i][j];

	for(int i=0; i<triang_num;i++)
		for(int j=0; j<3; j++)
			d.faces[i][j] = temp_f[d_Order[i]][j];
			    
    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Fish!");
    // glClearColor(1.0f,1.0f,1.0f,0.0f);      // Grey ackround color
	glMatrixMode(GL_PROJECTION);
	// glClear(GL_COLOR_BUFFER_BIT);

	// glEnable( GL_BLEND );
	// glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
    glutDisplayFunc(render);
    glutMainLoop();
	
	fclose(file);
    fclose(file2);
	fclose(file3);
    fclose(file4);
    return 0;
}