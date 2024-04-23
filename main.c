#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#define IX(i,j,k) ((i)+(M+2)*(j) + (M+2)*(N+2)*(k)) 
#define MAX(a,b)            (((a) > (b)) ? (a) : (b))

#define WINDOW_TITLE "Fluid"
#define WINDOW_WIDTH 768
#define WINDOW_HEIGHT 768
#define SIZE 42 // Best not to raise this very high

extern void dens_step ( int M, int N, int O, float * x, float * x0, float * u, float * v, float * w, float diff, float dt );
extern void vel_step (int M, int N, int O, float * u, float * v,  float * w, float * u0, float * v0, float * w0, float visc, float dt );


//fluid field information
static int M = SIZE; // grid x
static int N = SIZE; // grid y
static int O = SIZE; // grid z
static float dt = 0.4f; // time delta
static float diff = 0.0f; // diffuse
static float visc = 0.0f; // viscosity
static float force = 10.0f;  // added on keypress on an axis
static float source = 200.0f; // density
static float source_alpha =  0.05; //for displaying density

static int addforce[3] = {0, 0, 0};
static int addsource = 0;

static float * u, * v, *w, * u_prev, * v_prev, * w_prev;
static float * dens, * dens_prev;

static int dvel = 0;
static int dhelp = 1;
static int daxis = 1;

static int win_id;
static int win_x = WINDOW_WIDTH; 
static int  win_y = WINDOW_HEIGHT;
static int mouse_down[3];
static int mx, my;

GLfloat trans[3];
GLfloat azimuth = -45.f, elevation = 30.f;

static void free_data ( void )
{
	if ( u ) free ( u );
	if ( v ) free ( v );
	if ( w ) free ( w );
	if ( u_prev ) free ( u_prev );
	if ( v_prev ) free ( v_prev );
	if ( w_prev ) free ( w_prev );
	if ( dens ) free ( dens );
	if ( dens_prev ) free ( dens_prev );
}

static void clear_data ( void )
{
	int i, size=(M+2)*(N+2)*(O+2);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = w[i] = u_prev[i] = v_prev[i] =w_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}

	addforce[0] = addforce[1] = addforce[2] = 0;
}

static int allocate_data ( void )
{
	int size = (M+2)*(N+2)*(O+2);

	u			= (float *) malloc ( size*sizeof(float) );
	v			= (float *) malloc ( size*sizeof(float) );
	w			= (float *) malloc ( size*sizeof(float) );
	u_prev		= (float *) malloc ( size*sizeof(float) );
	v_prev		= (float *) malloc ( size*sizeof(float) );
	w_prev		= (float *) malloc ( size*sizeof(float) );
	dens		= (float *) malloc ( size*sizeof(float) );	
	dens_prev	= (float *) malloc ( size*sizeof(float) );

	if ( !u || !v || !w || !u_prev || !v_prev || !w_prev || !dens || !dens_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}

	return ( 1 );
}

static void get_force_source ( float * d, float * u, float * v, float * w )
{
	int i, j, k, size=(M+2)*(N+2)*(O+2);;

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = w[i]= d[i] = 0.0f;
	}

	if(addforce[0]==1) // x
	{
		i=2,
		j=N/2;
		k=O/2;

		if ( i<1 || i>M || j<1 || j>N || k <1 || k>O) return;
		u[IX(i,j,k)] = force*10;
		addforce[0] = 0;
	}	

	if(addforce[1]==1)
	{
		i=M/2,
		j=2;
		k=O/2;

		if ( i<1 || i>M || j<1 || j>N || k <1 || k>O) return;
		v[IX(i,j,k)] = force*10;
		addforce[1] = 0;
	}	

	if(addforce[2]==1) // y
	{
		i=M/2,
		j=N/2;
		k=2;

		if ( i<1 || i>M || j<1 || j>N || k <1 || k>O) return;
		w[IX(i,j,k)] = force*10; 	
		addforce[2] = 0;
	}	

	if(addsource==1)
	{
		i=M/2,
		j=N/2;
		k=O/2;
		d[IX(i,j,k)] = source;
		addsource = 0;
	}
	
	return;
}

static void draw_velocity ( void )
{
	int i, j, k;
	float x, y, z, h;

	h = 1.0f/MAX(MAX(M, N), MAX(N, O));

	glColor3f ( 1.0f, 1.0f, 1.0f );
	glLineWidth ( 1.0f );

	glBegin ( GL_LINES );

	for ( i=1; i<=M; i++ ) {
		x = (i-0.5f)*h;
		for ( j=1; j<=N; j++ ) {
			y = (j-0.5f)*h;
			for ( k=1; k<=O; k++ ) {
				z = (k-0.5f)*h;

				glVertex3f ( x, y, z );
				glVertex3f ( x+u[IX(i,j,k)], y+v[IX(i,j,k)], z+w[IX(i,j,k)] );
			}
		}
	}

	glEnd ();
}

static void draw_axis()
{

	glLineWidth ( 1.0f );
	glBegin (GL_LINES);	

	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f (0.25f, 0.f, 0.25f);
	glVertex3f (1.0f, 0.f, 0.25f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f (0.25f, 0.f, 0.25f);
	glVertex3f (0.25f, 1.0f, 0.25f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f (0.25f, 0.f, 0.25f);
	glVertex3f (0.25f, 0.f, 1.0f);

	glEnd();
}


static void draw_density ( void )
{
    int i, j, k;
	float x, y,z, h, d000, d010, d100, d110,d001, d011, d101, d111;
	
	h = 1.0f/MAX(MAX(M, N), MAX(N, O));

	glBegin ( GL_QUADS );

	for ( i=0; i<=M; i++ ) {
		x = (i-0.5f)*h;
		for ( j=0; j<=N; j++ ) {
			y = (j-0.5f)*h;
			for ( k=0; k<=O; k++ ) {
				z = (k-0.5f)*h;

				d000 = dens[IX(i,j,k)];
				d010 = dens[IX(i,j+1,k)];
				d100 = dens[IX(i+1,j,k)];
				d110 = dens[IX(i+1,j+1,k)];

				d001 = dens[IX(i,j,k+1)];
				d011 = dens[IX(i,j+1,k+1)];
				d101 = dens[IX(i+1,j,k+1)];
				d111 = dens[IX(i+1,j+1,k+1)];				
                
				// draw density as a cube of quads

				glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
				glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
				glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );
				glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );

				glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
				glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
				glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );
				glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );

				glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
				glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
				glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );
				glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
                
				glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
				glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
				glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
				glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );

				glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );
				glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
				glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );
				glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );

				glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
				glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
				glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
				glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h, y+h, z+h );				
			}
		}
	}

	glEnd ();
}

int init(void)
{
	if ( !allocate_data () ) 
		return (0);
	clear_data ();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0);

	return (1);
}

void draw_text(GLint x, GLint y, char* s, GLfloat r, GLfloat g, GLfloat b)
{
	int lines;
	char* p;

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	gluOrtho2D(0, win_x, 0, win_y);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glColor3f(r,g,b);
	glRasterPos2i(x, y);
	for(p = s, lines = 0; *p; p++) {
		if (*p == '\n')
		{
			lines++;
			glRasterPos2i(x, y-(lines*24));
		}
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *p);
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

void renderBitmapString(float x, float y ,const char *string){
    const char *c;
    glRasterPos2f(x, y);
    for (c=string; *c != '\0'; c++) {
        if (*c=='\n')
            glRasterPos2f(x, y-=0.12);
        else
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);
    }
}

void draw_help()
{
    glColor4f(0.8f, 0.8f, 0.8f, 1.0f);
    renderBitmapString(-0.5, 0.8,
                     "Help:\n" \
                     "      'X' key - add source at center\n" \
                     "      'W' key - apply force x-axis\n" \
                     "      'D' key - apply force y-axis\n" \
                     "      'S' key - apply force z-axis\n" \
                     "      'C' key - clear simulation\n" \
                     "      'V' key - show/hide velocities\n" \
                     "      'A' Key - show/hide the XYZ axis\n" \
                     "      'H' key - show/hide this help message\n"\
                     "      Left click  - pan from location\n" \
                     "      Right click - rotate cube\n"\
                     "      ESC - quit"\
                    );

}



void sim_main(void)
{

	get_force_source( dens_prev, u_prev, v_prev, w_prev );
	vel_step ( M,N,O, u, v, w, u_prev, v_prev,w_prev, visc, dt );
	dens_step ( M,N,O, dens, dens_prev, u, v, w, diff, dt );	

}

int shutdown(void)
{
	
	free_data ();
	
	return 0;
} 

void sim_reset()
{
	clear_data();
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glTranslatef(trans[0], trans[1], trans[2]);
	glRotatef(elevation, 1.0f, 0.0f, 0.0f);
	glRotatef(azimuth, 0.0f, 1.0f, 0.0f);

	// toggle display modes
	if ( dvel ) draw_velocity ();
	else		draw_density ();
	if (dhelp) draw_help();
	if (daxis) draw_axis();

	glEnd();
	glPopMatrix();
	glFlush();
	glutSwapBuffers();
}

static void reshape_func ( int width, int height )
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (float)width/height, 0.1, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);

	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}


static void key_func ( unsigned char key, int x, int y )
{
	
		switch (key) {
			case 27:		// ESC key
			    free_data ();
			    exit ( 0 );
				break;
			case 'w':       // 'W' key - apply force x-axis
				addforce[1] = 1;
				break;
			case 'd':       // 'D' key - apply force y-axis
				addforce[0] = 1;
				break;
			case 's':       // 'S' key - apply force z-axis
				addforce[2] = 1;
				break;
			case 'x':       // 'X' key - add source at centre
				addsource = 1;
				break;
			case 'c':       // 'C' key - clear simulation
				sim_reset();
				break;
			case 'v':       // 'V' key show velocity
				dvel = !dvel;  // toggle show velocity
				break;
			case 'a':       // 'A' Key draw axis
				daxis = !daxis;  // toggle draw axis
				break;
			case 'h':       // 'H' key - get help*/
				dhelp = !dhelp;  // toggle show help
				break;
		}
}

static void mouse_func( int button, int state, int x, int y ) {
	mx = x;
	my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y ) {
	// when move mouse
	float dx = x - mx;
	float dy = y - my;

	if( mouse_down[0] ) {
		azimuth += dx * 0.1f;
		elevation += dy * 0.1f;
	}

	if( mouse_down[1] ) {
		trans[0] -= dx / 1500.0f;
		trans[1] -= dy / 1500.0f;
	}
	
	if( mouse_down[2] ) {
		trans[2] -= (dx+dy) / 100.0f;
	}

	mx = x;
	my = y;
}

static void idle_func (void ) {
	sim_main();
	
	glutPostRedisplay();
	
}


static void open_glut_window ( void )
{
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( win_x, win_y );
	win_id = glutCreateWindow ( "3D Fluid Simulation" );

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	glutKeyboardFunc ( key_func );
	glutMouseFunc ( mouse_func );
	glutMotionFunc ( motion_func );
	glutReshapeFunc ( reshape_func );
	glutIdleFunc ( idle_func );
	glutDisplayFunc ( display);
}


int main ( int argc, char ** argv )
{

	glutInit ( &argc, argv );	
	open_glut_window();
	init();
	glutMainLoop();
	shutdown();

	return 0;
}
