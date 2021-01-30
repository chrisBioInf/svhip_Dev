#include <stdio.h>
#include <libtcc.h>


const char* tree = "\
#include <stdlib.h> \
#include <stdio.h> \
#include <math.h> \
int foo(float features[]) { \
    int classes[2]; \
    if (features[2] <= 0.15435723215341568) { \
        if (features[0] <= 0.21741315722465515) {\
            if (features[0] <= -0.1797432079911232) {\
                classes[0] = 0; \
                classes[1] = 125; \  
            } else { \
                if (features[1] <= -0.2515578716993332) {\
                    classes[0] = 9; \
                    classes[1] = 0; \
                } else {\
                    if (features[1] <= 0.1522204354405403) {\
                        if (features[0] <= -0.008525954326614738) {\
                            if (features[1] <= -0.0260246186517179) {\
                                classes[0] = 0; \
                                classes[1] = 3; \
                            } else {\
                                if (features[2] <= -0.39103659242391586) {\
                                    classes[0] = 0; \
                                    classes[1] = 1; \
                                } else {\
                                    classes[0] = 7;\ 
                                    classes[1] = 0; \
                                }\
                            }\
                        } else {\
                            if (features[2] <= 0.06833194009959698) {\
                                classes[0] = 0; \
                                classes[1] = 19;\ 
                            } else {\
                                classes[0] = 2;\ 
                                classes[1] = 0; \
                            }\
                        }\
                    } else {\
                        classes[0] = 0;\ 
                        classes[1] = 58;\ 
                    }\
                }\
            }\
        } else {\
            classes[0] = 30; \
            classes[1] = 0; \
        }\
    } else {\
        if (features[0] <= -0.41063331067562103) {\
            classes[0] = 0; \
            classes[1] = 1; \
        } else {\
            classes[0] = 145; \
            classes[1] = 0; \
        }\
    }\
    int class_idx = 0;\
    int class_val = classes[0];\
    int i;\
    for (i = 1; i < 2; i++) {\
        if (classes[i] > class_val) {\
            class_idx = i;\
            class_val = classes[i];\
        }\
    }\
    return class_idx;}" 
;


int main(int argc, char **argv)
{
    TCCState* s = tcc_new();
    if(!s){
        printf("Can’t create a TCC context\n");
        return 1;
    }
    tcc_set_output_type(s, TCC_OUTPUT_MEMORY);

    if (tcc_compile_string(s, tree) != 0) {
        printf("Compilation error !\n");
        return 2;
    }
    tcc_relocate(s, TCC_RELOCATE_AUTO);

    int (*const foo)(float features[]) = tcc_get_symbol(s, "foo");

    /* Read in some test values and predict them: */

    float vector1[3] = {-0.157, 0.31, -0.26};
    float vector2[3] = {-0.18, 0.56, -0.334};
    float vector3[3] = {-0.24, 0.08, -0.47};
    float vector4[3] = {0.37, 0.2, -0.81};
    float vector5[3] = {-0.25, 0.39, -0.41};

    int ret_val1 = foo(vector1);
    printf("Predicted class: %d \n", ret_val1);
    int ret_val2 = foo(vector2);
    printf("Predicted class: %d \n", ret_val2);
    int ret_val3 = foo(vector3);
    printf("Predicted class: %d \n", ret_val3);
    int ret_val4 = foo(vector4);
    printf("Predicted class: %d \n", ret_val4);
    int ret_val5 = foo(vector5);
    printf("Predicted class: %d \n", ret_val5);

    tcc_delete(s);
    return 0;
}
