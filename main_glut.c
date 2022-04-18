#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

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

void read_data( int M, int N, float m[M][N],  FILE *file){

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
	/*GLfloat v0[3] =  {-0.52, 0.34 -0.04};
	GLfloat v1[3] =  {-0.246, -0.65, -0.15};
	GLfloat v2[3] =  {-0.15, -0.50, 0.0};
	GLfloat c0[3] =  {0.742323, 0.229279, 0.356321};
	GLfloat c1[3] =  {0.356321, 0.364300, 0.403574};
	GLfloat c2[3] =  {0.403574, 0.321502, 0.000000};
        
	for(int i = 0; i<3; i++){
		v.v0[i] = v0[i];
		v.v1[i] = v1[i];
		v.v2[i] = v2[i];
		c.c0[i] = c0[i];	
		c.c1[i] = c1[i];
		c.c2[i] = c2[i];
    }*/	
	
	//c.c0[0] = (GLfloat)rand() / (GLfloat)RAND_MAX;
	//c.c0[1] = (GLfloat)rand() / (GLfloat)RAND_MAX;
	//c.c0[2] = (GLfloat)rand() / (GLfloat)RAND_MAX;
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
	
	float verts2d[4999][2];
	float vcolors[4999][3];
	float faces[10000][3];
	float depth[4999][1];

	read_data(4999,2, verts2d,file);
	read_data(4999, 3, vcolors,file2);
	read_data(10000, 3, faces, file3);
	read_data(4999, 1, depth, file4);
	
    

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