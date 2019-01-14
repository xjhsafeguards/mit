#ifndef INPUT_H
#define INPUT_H

#include "gheader.h"

class Input{
public:
  Input(){}

  void init();
  void iarg(int argc,char** argv);

};

extern Input *INPUT;

#endif
