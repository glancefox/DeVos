#include "DeVos.h"
#include <QtWidgets/QApplication>


int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	DeVos w;
	w.show();
	return a.exec();

}
