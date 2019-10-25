/*
 *  MRAG_IO_ArgumentParser.h
 *  MRAG
 *
 *	This argument parser assumes that all arguments are optional ie, each of
 *the argument names is preceded by a '-' all arguments are however NOT optional
 *to avoid a mess with default values and returned values when not found!
 *
 *	More converter could be required:
 *		add as needed
 *			TypeName as{TypeName}() in Value
 *
 *  Created by Christian Conti on 6/7/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#ifndef ARGUMENTPARSER_H_
#define ARGUMENTPARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <string>
#include <vector>

using namespace std;

namespace MRAG {

class Value {
 private:
  string content;

 public:
  Value() : content("") {}

  Value(string content_)
      : content(content_) { /*printf("%s\n",content.c_str());*/
  }

  double asDouble() const {
    if (content == "") return 0;
    return (double)atof(content.c_str());
  }

  int asInt() const {
    if (content == "") return 0;
    return atoi(content.c_str());
  }

  bool asBool() const {
    if (content == "") return false;
    if (content == "0") return false;
    if (content == "false") return false;

    return true;
  }

  string asString() const { return content; }
};

class ArgumentParser {
 private:
  map<string, Value> mapArguments;

  const int iArgC;
  const char **vArgV;
  bool bStrictMode;

 public:
  Value operator()(const string arg) {
    if (bStrictMode) {
      map<string, Value>::const_iterator it = mapArguments.find(arg);

      if (it == mapArguments.end()) {
        printf("Runtime option NOT SPECIFIED! ABORTING! name: %s\n",
               arg.data());
        abort();
      }
    }

    printf("%s is %s\n", arg.data(), mapArguments[arg].asString().data());
    return mapArguments[arg];
  }

  ArgumentParser(const int argc, const char **argv)
      : mapArguments(), iArgC(argc), vArgV(argv), bStrictMode(false) {
    for (int i = 1; i < argc; i++)
      if (argv[i][0] == '-') {
        string values = "";
        int itemCount = 0;

        for (int j = i + 1; j < argc; j++)
          if (argv[j][0] == '-')
            break;
          else {
            if (strcmp(values.c_str(), "")) values += ' ';

            values += argv[j];
            itemCount++;
          }

        if (itemCount == 0) values += '1';
        mapArguments[argv[i]] = Value(values);
        i += itemCount;
      }

    printf("found %ld arguments of %d\n", mapArguments.size(), argc);
  }

  int getargc() const { return iArgC; }

  const char **getargv() const { return vArgV; }

  void set_strict_mode() { bStrictMode = true; }

  void unset_strict_mode() { bStrictMode = false; }

  void save_options(string path = ".") {
    string options;
    for (map<string, Value>::const_iterator it = mapArguments.begin();
         it != mapArguments.end(); it++) {
      options += it->first + " " + it->second.asString() + " ";
    }
    string filepath = (path + "/" + string("argumentparser.log"));
    FILE *f = fopen(filepath.data(), "a");
    if (f == NULL) {
      printf("impossible to write %s.\n", filepath.data());
      return;
    }
    fprintf(f, "%s\n", options.data());
    fclose(f);
  }
};
}  // namespace MRAG

#endif
