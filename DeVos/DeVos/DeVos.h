#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_DeVos.h"

class DeVos : public QMainWindow
{
	Q_OBJECT

public:
	DeVos(QWidget *parent = Q_NULLPTR);

private:
	Ui::DeVosClass ui;
};
