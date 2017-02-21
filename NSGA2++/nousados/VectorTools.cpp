/*
 * VectorTools.cpp
 *
 *  Created on: 06-ago-2009
 *      Author: antonio
 */

#include <iostream>
#include <sstream>
#include "VectorTools.h"

using namespace std;

vector<vector<double> > VectorTools::loadDelimVector(string fileName, string delims)
{
   //regular expression:
   //  zero or more spaces followed by
   //  one or more delimiters followed by
   //  zero or more spaces.
   boost::regex reg("\\s*[" + delims + "]+\\s*");
   boost::regex rex_comment("\\s*[//|#].*");

   ifstream dataFile(fileName.c_str(), ios::in);
   if (dataFile.fail())
      throw std::exception();

   vector<vector<double> > matrix;
   string line;

   while ( getline(dataFile, line) )
   {
      if (line.size() == 0 || boost::regex_match(line, rex_comment))
         continue;

      vector<double> vec;
      const boost::sregex_token_iterator end;
      for (boost::sregex_token_iterator iter(line.begin(), line.end(), reg, -1); iter != end; ++iter)
      {
         try {
             if (string(*iter).length() > 0)
                vec.push_back( boost::lexical_cast<double>(string(*iter)) );
         }
         catch ( const boost::bad_lexical_cast &exc )
         {   // conversion failed, exception thrown by lexical_cast and caught
            cerr << "\ntoken: '" << string(*iter) << "'" << endl; 
            throw;
         }
      }

      matrix.push_back(vec);
   }

   dataFile.close();

   return matrix;
}

vector<valarray<double> > VectorTools::loadDelimArray(string fileName, string delims)
{
   //regular expression:
   //  zero or more spaces followed by
   //  one or more delimiters followed by
   //  zero or more spaces.
   boost::regex reg("\\s*[" + delims + "]+\\s*");
   boost::regex rex_comment("\\s*[//|#].*");

   ifstream dataFile(fileName.c_str(), ios::in);
   if (dataFile.fail())
      throw std::exception();

   vector<valarray<double> > matrix;
   string line;

   while ( getline(dataFile, line) )
   {
      if (line.size() == 0 || boost::regex_match(line, rex_comment))
         break;

      vector<double> vec;
      const boost::sregex_token_iterator end;
      for (boost::sregex_token_iterator iter(line.begin(), line.end(), reg, -1); iter != end; ++iter)
      {
         try {
            if (string(*iter).length() > 0)
               vec.push_back( boost::lexical_cast<double>(string(*iter)) );
         }
         catch ( const boost::bad_lexical_cast &exc )
         {   // conversion failed, exception thrown by lexical_cast and caught
            throw;
         }
      }

      //matrix.push_back(valarray<double>(vec.data(), vec.size()));
      matrix.emplace_back(vec.data(), vec.size());
   }

   dataFile.close();

   return matrix;
}
