#ifndef CUSTOM_H_INCLUDED
#define CUSTOM_H_INCLUDED

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* MATH STUFF */
/* computes the exponent of a double precision floating point number in scientific notation */
int exponent(double f){
    double temp;
    if((f != 0.0) || (f != -0.0)){
        temp = log10(fabs(f));
        if(temp < 0)
            temp -= 1.0;
        }
    else
        temp = 0.0;
    return (int) temp;
}

/* computes the mantissa (significant digits) of a double precision floating point number in scientific notation */
double mantissa(double f){
    double temp;
    if((f != 0.0) || (f != -0.0))
        temp = f*pow(10,-exponent(f));
    else
        temp = 0.0;
    return temp;
}

/* prints a number in Mathematica readable format to a file */
void fprintM(FILE* file_name, double f){
    double tempm;
    double tempe;
    if((f != 0.0) || (f != -0.0)){
        // computes the exponent
        tempe = log10(fabs(f));
        if(tempe < 0)
            tempe -= 1.0;
        // computes the mantissa
        tempm = f*pow(10,-(int)tempe);
    }
    else{
        tempm = 0.0;
        tempe = 0.0;
    }
    // prints output to file
    fprintf(file_name, "%lf*^%i", tempm, (int)tempe);
}

/* STRING MANIPULATION */

/* replaces all instances of orig[] in in[] with rep[] */
void replace_string(char in[], char orig[], char rep[]){
    char temp[1000] = {};    // not sure how to make this dynamic, hopefully this is large enough
    int len_in = strlen(in), len_orig = strlen(orig), len_rep = strlen(rep);    // lengths of the strings involved
    int i, j, k, l, m=0;    // i,j,k are iterators, l counts character matching, m counts number of replacements

    // iterates over the input string
    for(i=0, k=0; i<len_in; i++, k++){
        // triggers when character in input matches first character of the string to be replaced
        if(in[i] == orig[0] && i <= (len_in - len_orig)){
            l = 0;
            // checks to see if each successive character matches
            for(j=0; j<len_orig; j++){
                if(in[i+j] == orig[j])
                    l++;
                else
                    break;
            }
            // triggers if a string match is found
            if(l == len_orig){
                // adds the replacement string to the temporary string
                for(j=0; j<len_rep; j++, k++){
                    temp[k] = rep[j];
                }
                k--;    // need to move this back so it can be iterated correctly by the enclosing for loop
                i += len_orig - 1;
                m++;    // counts the number of times the replacement has been made
            }
            else{
                temp[k] = in[i];
            }
        }
        else{
            temp[k] = in[i];
        }
    }
    for(i=0; i<len_in + m*(len_rep - len_orig); i++){
        strcpy(in,temp);    // overwrites the old string with the new one
    }
}


#endif // CUSTOM_H_INCLUDED
