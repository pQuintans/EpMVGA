#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <list>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <time.h>
//#include <functional>

#include "of.h"
#include "ofOffPointsReader.h"
#include "Handler.hpp" 
#include "GL_Interactor.h"
#include "ColorRGBA.hpp"
#include "Cores.h"
#include "Point.hpp"
#include "printof.hpp"


#include "CommandComponent.hpp"
#include "MyCommands.hpp"

#include "ofVertexStarIteratorSurfaceVertex.h"


clock_t start_insert;
clock_t end_insert;
clock_t start_print;
clock_t end_print;



using namespace std;
using namespace of;

//Define o tamanho da tela.
scrInteractor *Interactor = new scrInteractor(900, 650);

//Define a malha a ser usada.
typedef of::MyofDefault2D TTraits;
typedef of::ofMesh<TTraits> TMesh;
TMesh *malha;
Handler<TMesh> meshHandler;

typedef PrintOf<TTraits> TPrintOf;

TPrintOf *Print;

typedef MyCommands<TPrintOf> TMyCommands;
typedef CommandComponent TAllCommands;

ofVtkWriter<TTraits> writer;
TAllCommands *allCommands;

//##################################################################//

double* p = new double[2];
double initialCell = 167;
int timesExecuted = 0;

//##################################################################//

double* SubtractCoordinates(double* a, double* b);
double CrossMutiplyAndDiferenceCoordinate(double* a, double* b);
int CalculateDirection(of::ofMyCell<of::MyofDefault2D>* cell);

////////////////////////////////////////////////////////////////////////
int type = 3;
//CASO 1 EXECUTA CRUST
//CASO 2 EXECUTA BETA-SKELETON
//CASO 3 EXECUTA ARVORE
////////////////////////////////////////////////////////////////////////

double* SubtractCoordinates(double* a, double* b) {
    double* coordinate = new double[2];
    coordinate[0] = a[0] - b[0];
    coordinate[1] = a[1] - b[1];
    return coordinate;
}

double CrossMutiplyAndDiferenceCoordinate(double* a, double* b) {
    return (a[1] * b[0]) - (a[0] * b[1]);
}

int CalculateDirection(of::ofMyCell<of::MyofDefault2D>* cell){
    double* a = malha->getVertex(cell->getVertexId(0))->getCoords();
    double* b = malha->getVertex(cell->getVertexId(1))->getCoords();
    double* c = malha->getVertex(cell->getVertexId(2))->getCoords();


    double* diferenceBA = SubtractCoordinates(b, a);
    double* diferencePA = SubtractCoordinates(p, a);
    double* diferenceCA = SubtractCoordinates(c, a);

    double gama =  CrossMutiplyAndDiferenceCoordinate(diferencePA, diferenceBA) / CrossMutiplyAndDiferenceCoordinate(diferenceCA, diferenceBA);
    double beta = (diferencePA[0] - gama * diferenceCA[0]) / diferenceBA[0];
    double alpha = 1 - beta - gama;

    if(alpha > 0 && beta > 0 && gama > 0)
        return -2;

    if(alpha < beta && alpha < gama) {
        return cell->getMateVertexId(cell->getVertexId(0));
    } else if(beta < alpha && beta < gama) {
        return cell->getMateVertexId(cell->getVertexId(1));
    } else {
        return cell->getMateVertexId(cell->getVertexId(2));
    }

    return -2;
}

void RenderScene(void){
	allCommands->Execute();
    Print->Vertices(malha,blue,3);

    of::ofMyCell<of::MyofDefault2D>* cell;
    int index = initialCell;
    int lastIndex;

    Print->FacesWireframe(malha,grey,3);
    Print->Face(malha->getCell(initialCell), green);

    if(Interactor->mouseLeftIsClicked() == true){
        timesExecuted++;
        p[0] = Interactor->getP()[0];
        p[1] = Interactor->getP()[1];

        while(index > -1) {
            cell = malha->getCell(index);
            lastIndex = index;
            index = CalculateDirection(cell);

            if(index > -1)
                Print->Face(malha->getCell(lastIndex), grey);
        }

        if(index == -1) {
            std::cout<< std::endl<< "(" << timesExecuted << ") Nao e possivel acessar essa celula" <<std::endl<<std::endl;
        } else {
            std::cout<< std::endl<< "(" << timesExecuted << ") Celula encontrada" <<std::endl<<std::endl;
        }

        Print->Face(malha->getCell(lastIndex), (index == -2 ? red : grey));

        glFinish();
        glutSwapBuffers();
    }

	glFinish();
	glutSwapBuffers();
}



void HandleKeyboard(unsigned char key, int x, int y){
	double coords[3];
	char *xs[10];
	allCommands->Keyboard(key);
	
	switch (key) {

		case 'e':
			exit(1);
		break;
		case 'v':
			coords[0]=x;
			coords[1]=-y;
			coords[2]=0.0;
			malha->addVertex(coords);
		break;
		case 's':
			
			
		break;

		case 'd':
			
			
		break;
	

	}
    
	
    Interactor->Refresh_List();
	glutPostRedisplay();

}

using namespace std;

int main(int *argc, char **argv)
{

  ofRuppert2D<MyofDefault2D> ruppert;
  ofPoints2DReader<MyofDefault2D> reader;
  ofVtkWriter<MyofDefault2D> writer;
  Interactor->setDraw(RenderScene);
	meshHandler.Set(new TMesh());
      char *fileBrasil = "../Brasil.off";

     
    reader.readOffFile(fileBrasil);
    
    ruppert.execute2D(reader.getLv(),reader.getLids(),true);
    //writer.write(ruppert.getMesh(),"out.vtk",reader.getNorma(),ruppert.getNumberOfInsertedVertices());
  
  meshHandler = ruppert.getMesh();
  malha = ruppert.getMesh();
  
  
  Print = new TPrintOf(meshHandler);

	allCommands = new TMyCommands(Print, Interactor);

	double a,x1,x2,y1,y2,z1,z2; 

	of::ofVerticesIterator<TTraits> iv(&meshHandler);

	iv.initialize();
	x1 = x2 = iv->getCoord(0);
	y1 = y2 = iv->getCoord(1);
	z1 = z2 = iv->getCoord(2);

	for(iv.initialize(); iv.notFinish(); ++iv){
		if(iv->getCoord(0) < x1) x1 = a = iv->getCoord(0);
		if(iv->getCoord(0) > x2) x2 = a = iv->getCoord(0);
		if(iv->getCoord(1) < y1) y1 = a = iv->getCoord(1);
		if(iv->getCoord(1) > y2) y2 = a = iv->getCoord(1);
		if(iv->getCoord(2) < z1) z1 = a = iv->getCoord(2);
		if(iv->getCoord(2) > z2) z2 = a = iv->getCoord(2);
	}

	double maxdim;
	maxdim = fabs(x2 - x1);
	if(maxdim < fabs(y2 - y1)) maxdim = fabs(y2 - y1);
	if(maxdim < fabs(z2 - z1)) maxdim = fabs(z2 - z1);

	maxdim *= 0.6;
	
	Point center((x1+x2)/2.0, (y1+y2)/2.0, (y1+y2)/2.0 );
	Interactor->Init(center[0]-maxdim, center[0]+maxdim,
					center[1]-maxdim, center[1]+maxdim,
					center[2]-maxdim, center[2]+maxdim,argc,argv);

	
	
    AddKeyboard(HandleKeyboard);

	allCommands->Help(std::cout);
    std::cout<< std::endl<< "Press \"?\" key for help" <<std::endl<<std::endl;
	double t;
	
	Init_Interactor();

  
  return EXIT_SUCCESS;
}
