#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

//float canvas[512][512];

#define ALPHA 0.5
typedef struct 
{
	GLfloat v0[3];
	GLfloat v1[3];
	GLfloat v2[3];
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
	int x;
	int y;
	float rgb[3];
}point;


float interpolate_color(point p1, point p2, point p){
	float r1  = sqrt(pow(p2.x - p.x,2) + pow(p2.y - p.y,2));
	float r2  = sqrt(pow(p.x - p1.x,2) + pow(p.y - p1.y,2));
	float percent = r1 / (r1+r2);
	for(int i=0; i<3; i++){
		p.rgb[i] = percent * p1.rgb[i] + (1 - percent) * p2.rgb[i];
	}
}


void read_data_float( int M, int N, float m[M][N],  FILE *file){

	char line[50];

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

	char line[50];

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
		//printf("n = %d \n",n);
		// if(l=3)
		// 	printf("a[%d] = %d\n",l,o[l]);
		// 	printf("a[%d] = %d\n",r,o[r]);
		int p = partition2(a,o, n);
		q_sort(a,o,p);
		q_sort(&a[p+1],&o[p+1],n-p-1);
	}
}

void shade_triangle(triangle *v, color *c){
	
	/*for(int j =0; j<3; j++){
            v->v0[j] -=0.3;
            v->v1[j] +=0.3;
            v->v2[j] -=0.3;
    }*/
	GLfloat a = (GLfloat)rand() / (GLfloat)RAND_MAX;
	glBegin(GL_TRIANGLES);
		glColor4f(c->c0[0],c->c0[1],c->c0[2],a);
        glVertex3f(v->v1[0],v->v1[1],v->v1[2]);
		glColor4f(c->c1[0],c->c1[1],c->c1[2],a);
        glVertex3f(v->v1[0],v->v2[1],v->v2[2]);
		glColor4f(c->c2[0],c->c2[1],c->c2[2],a);
        glVertex3f(v->v2[0],v->v2[1],v->v2[2]);
    glEnd();

}

void render(void)
{
	triangle v;
	color c;
	glEnable( GL_BLEND );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
    glClearColor(0.5,0.5,0.5,0.3);      // Grey ackround color
    glClear(GL_COLOR_BUFFER_BIT);

	int facesLength = sizeof(d.faces) / sizeof(d.faces[0]);
	float *d_Median = (float *) malloc(facesLength * sizeof(float));
	for (int l = 0; l < facesLength; l++) d_Median[l] = 0.0;
	int *d_Order = (int *) malloc(facesLength * sizeof(int)); 
	for (int l = 0; l < facesLength; l++) d_Order[l] = l;

	for(int i=0; i<facesLength;i++){
		for(int j=0; j<3; j++){
			d_Median[i] = d.depth[d.faces[i][j]][0] + d_Median[i]; 
		}
		d_Median[i] = d_Median[i]/3;
	}

	float max = d_Median[0];
    float min = d_Median[0];

   
	
	q_sort(d_Median,d_Order,facesLength) ;
	
	for(int t=0; t<30; t++){
		for(int i = 0; i<3; i++){
			v.v0[i] = -1 + (2*((GLfloat)rand()) / (GLfloat)RAND_MAX);
			v.v1[i] = -1 + (2*((GLfloat)rand()) / (GLfloat)RAND_MAX);
			v.v2[i] = -1 + (2*((GLfloat)rand()) / (GLfloat)RAND_MAX);
			c.c0[i] = (GLfloat)rand() / (GLfloat)RAND_MAX;
			c.c1[i] = (GLfloat)rand() / (GLfloat)RAND_MAX;
			c.c2[i] = (GLfloat)rand() / (GLfloat)RAND_MAX;
			
		}		
		shade_triangle(&v, &c);
	}
glFlush();
glutPostRedisplay();
sleep(1);
    
}

int main(int argc, char** argv)
{	
	char const* const fileName = "verts2d.txt"; /* should check that argc > 1 */
	char const* const fileName2 = "vcolors.txt";
	char const* const fileName3 = "faces.txt";
	char const* const fileName4 = "depth.txt";

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

	read_data_int(4999,2, d.verts2d,file);
	read_data_float(4999, 3, d.vcolors,file2);
	read_data_int(10000, 3, d.faces, file3);
	read_data_float(4999, 1, d.depth, file4);
	//return 0;
	// printf("%d %d\n", d.verts2d[0][0], d.verts2d[0][1]);
	// printf("%d %d\n", d.verts2d[1][0], d.verts2d[1][1]);

	// printf("%d %d %d\n", d.faces[0][0], d.faces[0][1], d.faces[0][2]);
	// printf("%d %d %d\n", d.faces[1][0], d.faces[1][1], d.faces[1][2]);
	
    double rgb[3][3];
    float V1[3][2] = {{185.0, 272.0},
                      {355.0, 250.0},
                      {450.0, 300.0}};
    srand(time(0));

    for(int i = 0; i<3; i++){
        for(int j =0; j<3; j++){
            rgb[i][j] = (double)rand() / (double)RAND_MAX;
            //printf(" %f ", rgb[i][j]);
        }
        //printf("\n");
    }
       


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Fish!");
	
    //init();
    glutDisplayFunc(render);
    glutMainLoop();
	
	fclose(file);
    fclose(file2);
	fclose(file3);
    fclose(file4);
    return 0;
}