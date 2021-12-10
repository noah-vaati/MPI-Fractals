#ifndef _JULIA_H
#define _JULIA_H

/* Control Constants - feel free to play around with these settings to generate different images if you like. */

// Dimensions of the output image
#define WIDTH 800
#define HEIGHT 600
// Controls how close we "zoom in" on the fractal pattern (lower = closer)
#define SCALE 4.0
// Complex number contant that controls the pattern that's generated (try changing it!)
// Note that values closer to 0.0 may significantly increase the computation time.
// You can set either of the terms to anything in the range -1.0 to 1.0. Just make sure to
// keep the trailing 'I'.
#define C (-0.839 + 0.4I )
// Name of the file to same the output image in
#define FNAME "out.bmp"

/* Image Generation Constants - messing with these ones might break stuff */

#define DEPTH 24
#define MAX_COLOUR ((1L << DEPTH) - 1)
#define MAX_ITER ((1L << 16) - 1)
#define RES_FACTOR (sqrt(pow(WIDTH / 2, 2) + pow(HEIGHT / 2, 2)))
#define LIMIT 50
#define PI 3.14159265358979323846

/* Prototype Declarations */

void toRGB(unsigned int val,
           unsigned char *r, unsigned char *g, unsigned char *b);
unsigned int sum_array(unsigned int *array, int len);
void hist_eq(unsigned int *data, unsigned int *hist);
void write_bmp(unsigned int *data, char *fname);
unsigned int julia_iters(float complex z);
void compute_row(int row, unsigned int *data, unsigned int *hist);

#endif
