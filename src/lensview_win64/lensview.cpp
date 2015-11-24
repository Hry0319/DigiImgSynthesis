#include <cstdarg>
#include "lensdata.h"
#include "lensview.h"

// Global lenSystem object
static lensSystem lsystem;
// Global interface variables
static int lvWidth    = LVW;
static int lvHeight   = LVH;
static float lvZoom   = 1.0f;
static float lvTransX = 0.0f;
static float lvTransY = 0.0f;
static bool lvPressed[3] = {false,false,false};
static int lastx = 0, lasty = 0;

// Printf style formatting for rendering text in OpenGL
static void lvText(int x, int y, void *font, int height, 
				   char *format, ...) {
    va_list args; char buffer[256];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
    int w = lvWidth, h = lvHeight;

    glMatrixMode(GL_PROJECTION); 
	glPushMatrix(); glLoadIdentity();
    gluOrtho2D(0, w - 1, 0, h - 1);
    glMatrixMode(GL_MODELVIEW); 
	glPushMatrix(); glLoadIdentity();

    if ((x >= 0) && (y >= 0)) 
		glRasterPos2i(x,y);
    for (char *p = buffer; *p; ++p) 
		glutBitmapCharacter(font, *p);
   
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);  glPopMatrix();
}

static void applyTransforms(void) {
	glMatrixMode(GL_MODELVIEW);	 glLoadIdentity();
	glMatrixMode(GL_PROJECTION); glLoadIdentity();
	float orthoX = lvWidth / lvZoom;
	float orthoY = lvHeight / lvZoom;

	gluOrtho2D(orthoX, -orthoX, -orthoY, orthoY);
    glTranslated(lvTransX, lvTransY, 0.0f);
}

static void lvDisplay(void) {
    glClear(GL_COLOR_BUFFER_BIT);
	applyTransforms();

	glColor3f(0.25f, 0.25f, 0.25f);
	glLineWidth(2);
	drawLine(-lvWidth-lvTransX, 0, lvWidth-lvTransX, 0);
	drawLine(0, -lvHeight-lvTransY, 0, lvHeight-lvTransY);
	glLineWidth(1);

	lsystem.Draw();
	
	glMatrixMode(GL_PROJECTION); glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, lvWidth, lvHeight, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    drawSquare(0.25f, 0.25f, 0.25f, 0, 0, lvWidth, BARH);
	drawSquare(0.25f, 0.25f, 0.25f, 0, lvHeight, lvWidth, lvHeight-BARH);
	glLoadIdentity();

	gluOrtho2D(-lvWidth, lvWidth, lvHeight, -lvHeight);
	glColor3f(0, 1, 1);
	drawText(-lvWidth + 15, -lvTransY*lvZoom - 8, SBFONT, "+z");
	drawText(lvWidth - 50, -lvTransY*lvZoom - 8, SBFONT, "-z");
	drawText(-lvTransX*lvZoom, -lvHeight + 75, SBFONT, "+y");
	drawText(-lvTransX*lvZoom, lvHeight - 65, SBFONT, "-y");
	
    glColor3f(1, 1, 1);
	lvText(9, lvHeight-16, SBFONT, 0, "fstop %1.1f (aperture %1.1fmm) | "
		                              "object %1.1fmm | image %1.1fmm | mode: none", 
		                              0.,0.,0.,0.); 

	lvText(9, 8, SBFONT, 0, "F %1.1fmm | F' %1.1fmm | P %1.1fmm | "
		                    "P' %1.1fmm | exit pupil (%1.1fmm, %1.1fmm)", 
							0.,0.,0.,0.,0.,0.);
	glPopMatrix();
	glutSwapBuffers();
}

static void lvClick(int button, int state, int x, int y) {
	if (state == GLUT_DOWN)	
		lvPressed[button] = true;
	else if (state == GLUT_UP) 
		lvPressed[button] = false;

	lastx = x; lasty = y; 
	glutPostRedisplay();
}

static void lvMotion(int x, int y) {
	// Here we handle translation
	float norm = 1.0f / lvZoom;
	if (lvPressed[GLUT_LEFT_BUTTON]) {
		lvTransX += float(lastx - x)*2.0f*norm; 
        lvTransY += float(lasty - y)*2.0f*norm;
	}
	if (lvPressed[GLUT_MIDDLE_BUTTON]) {
		// Use this to change camera focus
	}
	// Here we handle zooming
	if (lvPressed[GLUT_RIGHT_BUTTON]) {
	    float zoom = 3.0f * float(x - lastx) / lvWidth;
	    if ((lvZoom + zoom) >= 1.0f)
			lvZoom += zoom;
	}

	lastx = x; lasty = y; 
	glutPostRedisplay();
}

static void lvResize(int w, int h) {
    lvWidth = w; lvHeight = h;
	glViewport(0, 0, lvWidth, lvHeight);
}

static void lvKey(unsigned char key, int x, int y) {
    switch (key) {
        case 27: /* Escape key */ exit(0);
		case 'm': // toggle trace mode
			break;
		case '-': // decrease fstop
			break;
		case '=': // increase fstop
			break;	
		default: break;
	}
	glutPostRedisplay();
}

int main(int argc, char** argv) {
	if (argc < 2) {
		printf("Usage: lensview <lens file>\n");
		exit (-1);
	}

	if (!lsystem.Load(argv[1])) {
		printf("Cannot open the lens file: %s.\n", argv[1]);
		exit(-1);
	}
	lvZoom = lvHeight / lsystem.maxAperture() / 3;
	lvZoom = maxv(1, lvZoom);

	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(lvWidth, lvHeight);
    glutCreateWindow("lensview");
    glutDisplayFunc(lvDisplay);
    glutMotionFunc(lvMotion);
    glutMouseFunc(lvClick);
    glutKeyboardFunc(lvKey);
	glutReshapeFunc(lvResize);
    glutIdleFunc(NULL);

	glMatrixMode(GL_MODELVIEW);  glLoadIdentity();
    glMatrixMode(GL_PROJECTION); glLoadIdentity();	
	glClearColor(BGR, BGG, BGB, 1);
	
	glutMainLoop();
	return (0);
}

	
	  
	  	  

