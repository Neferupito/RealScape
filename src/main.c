
#include <stdlib.h>
#include <stdio.h>
#include <stdalign.h> // Required for alignof in C11
#include <time.h>
#include <math.h>


int is_within_bounds(int i, int j, int rows, int cols) {
    return (i >= 0 && i< rows && j >= 0 && j < cols);
}

void write_xyz_to_binary_file(int rows, int cols, float x[rows][cols], float y[rows][cols], float z[rows][cols], char filename[13]) {
    // Open the file for binary writing
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Unable to open file for writing");
        return;
    }

    // Loop through the grid and write X, Y, Z values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float X = x[i][j];    // X coordinate
            float Y = y[i][j];    // Y coordinate
            float Z = z[i][j];       // Z value as a function of X and Y (you can modify this)

            // Write X, Y, Z as binary floats (4 bytes each)
            if (fwrite(&X, sizeof(float), 1, file) != 1 || 
                fwrite(&Y, sizeof(float), 1, file) != 1 || 
                fwrite(&Z, sizeof(float), 1, file) != 1) {
                perror("Error writing data to file");
                fclose(file);
                return;
            }
        }
    }

    // Close the file
    fclose(file);
    printf("Data written successfully to %s\n", filename);
}

void update_topo(int r, int c, float dx, float power_coef, float erosion_strength, float deposition_strength, float diffusion_coef, float z[r][c], float slope_map[r][c], float h_water[r][c]) {
  int i,j;
  for (i = 1; i < r-1; i++) {
      for (j = 1; j < c-1; j++) {
          z[i][j] += ((z[i][j+1] - 2*z[i][j] + z[i][j-1]) + (z[i+1][j] - 2*z[i][j] + z[i-1][j])) * diffusion_coef / pow(dx,2);
          if (slope_map[i][j] <= 0 && h_water[i][j] > 0) {
              z[i][j] += deposition_strength ;
          }
          else if (slope_map[i][j] > 0 && h_water[i][j] > 0) {
            if (h_water[i][j] > 1){h_water[i][j] = 1;}
              z[i][j] -=  pow(slope_map[i][j],power_coef)* pow(h_water[i][j],1)*erosion_strength;
          }
     }
 } 
} 


void Dinf_maps(int rows, int cols, float dx, float z[rows][cols], float dinf_slope_map[rows][cols], float dinf_split1[rows][cols], float dinf_split2[rows][cols], int dinf_split1_i[rows][cols], int dinf_split1_j[rows][cols],  int dinf_split2_i[rows][cols], int dinf_split2_j[rows][cols]) {
  const float PI = 3.1415926;
  float max_slope, r_max, e0, e1, e2, s1, s2, s, rg, r, alpha_high, alpha_low;
  int dinf_directions[9][2] = {{0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}};
  int dinf_directions_1[8][2] = {{0, 1}, {-1, 0}, {-1, 0}, {0, -1}, {0, -1}, {1, 0}, {1, 0}, {0, 1}};
  int dinf_directions_2[8][2] = {{-1, 1}, {-1, 1}, {-1, -1}, {-1, -1}, {1, -1}, {1, -1}, {1, 1}, {1, 1}};
  int ac[8] = {0, 1, 1, 2, 2, 3, 3, 4}, af[8] = {1, -1, 1, -1, 1, -1, 1, -1};
  int neighbour_id, i,j,k, k_max, i_1, j_1, i_2, j_2;

  // Update the grid
  for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
          dinf_slope_map[i][j] = -1;
          e0 = z[i][j];
          max_slope = 0;
          k_max = 0;
          r_max = 0;
          for (k = 0; k < 8; k++) {
              i_1 = i + dinf_directions_1[k][0];
              i_2 = i + dinf_directions_2[k][0];
              j_1 = j + dinf_directions_1[k][1];
              j_2 = j + dinf_directions_2[k][1];
              // Ensure the neighbor is within bounds
              if (is_within_bounds(i_1, j_1, rows, cols) && is_within_bounds(i_2, j_2, rows, cols)) {
                  e1 = z[i_1][j_1];
                  e2 = z[i_2][j_2];
                  s1 = (e0 - e1) / dx;
                  s2 = (e1 - e2) / dx;

                  r = atan2(s2 , s1);  // Calculate angle based on slope
                  s = sqrt(s1*s1 + s2*s2);  // Total slope

                  // Boundary handling for steepest slope limits
                  if (r < 0) {
                    r = 0; 
                    s = s1;
                    }
                  else if (r > atan2(dx , dx)) {
                    r = atan2(dx, dx); 
                    s = (e0 - e2) / (dx * sqrt(2));
                    }

                  if (s > max_slope) {
                    max_slope = s; 
                    k_max = k; 
                    r_max = r;
                  }
              }
          }
          // Flow redistribution
          if (max_slope > 0) {
              rg = af[k_max] * r_max + ac[k_max]* PI / 2 ;
              neighbour_id = rg / (PI/4);
              alpha_high = (PI/4)*(neighbour_id+1) - rg;
              alpha_low = rg - (PI/4)*neighbour_id;
              dinf_slope_map[i][j] = max_slope;
              dinf_split1[i][j] = alpha_low / (alpha_high + alpha_low);
              dinf_split2[i][j] = alpha_high / (alpha_high + alpha_low);
              dinf_split1_i[i][j] = i + dinf_directions[neighbour_id+1][0];
              dinf_split2_i[i][j] = i + dinf_directions[neighbour_id][0];
              dinf_split1_j[i][j] = j + dinf_directions[neighbour_id+1][1];
              dinf_split2_j[i][j] = j + dinf_directions[neighbour_id][1];

          }
      }
  }
  
}



void Dinf_update_flow(int rows, int cols, float rain, float h_water[rows][cols], float dinf_slope_map[rows][cols], float dinf_split1[rows][cols], float dinf_split2[rows][cols], int dinf_split1_i[rows][cols], int dinf_split1_j[rows][cols], int dinf_split2_i[rows][cols], int dinf_split2_j[rows][cols]) {
  int i,j;
  float (*h_water_tampon)[cols] = malloc(rows * sizeof *h_water_tampon);
  
  for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
          h_water_tampon[i][j] = h_water[i][j] + rain;
          h_water[i][j] = 0;
      }
  }

  for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
        if (dinf_slope_map[i][j] > 1.e-7) {
          h_water[dinf_split1_i[i][j]][dinf_split1_j[i][j]] += dinf_split1[i][j] * h_water_tampon[i][j];
          h_water[dinf_split2_i[i][j]][dinf_split2_j[i][j]] += dinf_split2[i][j] * h_water_tampon[i][j];
        }
      }
  }
  free(h_water_tampon);
}

void Dinf_step(int rows, int cols, float dx, float rain, float z[rows][cols], float h_water[rows][cols], float power_coef, float erosion_strength, float deposition_strength, float diffusion_coef) {
  float (*dinf_slope_map)[cols] = malloc(rows * sizeof *dinf_slope_map);
  float (*dinf_split1)[cols] = malloc(rows * sizeof *dinf_split1);
  float (*dinf_split2)[cols] = malloc(rows * sizeof *dinf_split2);
  int (*dinf_split1_i)[cols] = malloc(rows * sizeof *dinf_split1_i);
  int (*dinf_split1_j)[cols] = malloc(rows * sizeof *dinf_split1_j);
  int (*dinf_split2_i)[cols] = malloc(rows * sizeof *dinf_split2_i);
  int (*dinf_split2_j)[cols] = malloc(rows * sizeof *dinf_split2_j);
  
  Dinf_maps(rows, cols, dx, z, dinf_slope_map, dinf_split1, dinf_split2, dinf_split1_i, dinf_split1_j, dinf_split2_i, dinf_split2_j);
  Dinf_update_flow(rows, cols, rain, h_water, dinf_slope_map,dinf_split1, dinf_split2, dinf_split1_i, dinf_split1_j, dinf_split2_i, dinf_split2_j);
  update_topo(rows, cols, dx, power_coef, erosion_strength, deposition_strength, diffusion_coef, z, dinf_slope_map,  h_water);

  free(dinf_slope_map);
  free(dinf_split1);
  free(dinf_split2);
  free(dinf_split1_i);
  free(dinf_split1_j);
  free(dinf_split2_i);
  free(dinf_split2_j);
}

void D8_maps(int r, int c, float dx, float z[r][c], float d8_slope_map[r][c], int d8_i[r][c], int d8_j[r][c]) {
  int i,j;
  int neighbour_i, neighbour_j, neighbour_id;
  int max_i, max_j;
  float max_slope, neighbour_slope;
  float d8_distances[8] = {dx, dx * sqrt(2), dx, dx * sqrt(2), dx, dx * sqrt(2), dx, dx * sqrt(2)};
  int d8_directions[8][2] = {{-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};

  // Update the grid
  for (i = 0; i < r; i++) {
      for (j = 0; j < c; j++) {
          d8_slope_map[i][j] = -1;
          // Check the 8 neighbors
          max_slope = 0;
          for (neighbour_id = 0; neighbour_id < 8; neighbour_id++) {
              neighbour_i = i + d8_directions[neighbour_id][0];
              neighbour_j = j + d8_directions[neighbour_id][1];
              
              // Ensure the neighbor is within bounds
              if (is_within_bounds(neighbour_i, neighbour_j, r, c)) {
                  // Calculate the difference in elevation
                  neighbour_slope = (z[i][j] - z[neighbour_i][neighbour_j])/d8_distances[neighbour_id];
                  
                  // Find the steepest slope
                  if (neighbour_slope > max_slope) {
                      max_slope = neighbour_slope;
                      max_i = neighbour_i;
                      max_j = neighbour_j;
                  }
              }
          }
          if (max_slope > 0) {
              d8_i[i][j] = max_i;
              d8_j[i][j] = max_j;
              d8_slope_map[i][j] = max_slope;
          }
      }
  }
}

void D8_update_flow(int rows, int cols, float rain, float h_water[rows][cols], float d8_slope_map[rows][cols], int d8_i[rows][cols], int d8_j[rows][cols]) {
  int i,j;
  float (*h_water_tampon)[cols] = malloc(rows * sizeof *h_water_tampon);
  for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
          h_water_tampon[i][j] = h_water[i][j] + rain;
          h_water[i][j] = 0;
      }
  }

  for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
        if (d8_slope_map[i][j] > 0) {
          h_water[d8_i[i][j]][d8_j[i][j]] +=  h_water_tampon[i][j];
        }
      }
  }
  free(h_water_tampon);
}

void D8_step(int rows, int cols, float dx, float rain, float z[rows][cols], float h_water[rows][cols], float power_coef, float erosion_strength, float deposition_strength, float diffusion_coef) {
  float (*d8_slope_map)[cols] = malloc(rows * sizeof *d8_slope_map);
  int (*d8_i)[cols] = malloc(rows * sizeof *d8_i);
  int (*d8_j)[cols] = malloc(rows * sizeof *d8_j);
  
  D8_maps(rows, cols, dx, z, d8_slope_map, d8_i, d8_j);
  D8_update_flow(rows, cols, rain, h_water, d8_slope_map, d8_i, d8_j);
  update_topo(rows, cols, dx, power_coef, erosion_strength, deposition_strength, diffusion_coef, z, d8_slope_map,  h_water);

  free(d8_slope_map);
  free(d8_i);
  free(d8_j);
}
    
int main(int r, int c, int iter_max, float rain, float power_coef, float erosion_strength, float deposition_strength, float diffusion_coef) {
    
    // InOut
    char file_in[] = "../inout/in.xyz", file_out[] = "../inout/out.bin"; // Change this to your file name
    
    // Topo grid parameters
    int i, j, iter; // Loop counters
    int line; // Read file line number
    float x_file, y_file, z_file; // Store fscanf values
    float max_x, max_y, max_z, min_x, min_y, min_z; // Store min/max values
    float dx, dy; // Store dx and dy

    float (*x)[c] = malloc(r * sizeof *x);
    float (*y)[c] = malloc(r * sizeof *y);
    float (*z)[c] = malloc(r * sizeof *z);
    float (*h_water)[c] = malloc(r * sizeof *h_water);

    // Start the clock
    clock_t start_time = clock();

    // Set file_in pointer
    FILE *file_in_id;

    // Open a file in read mode
    file_in_id = fopen(file_in, "r");

    // Print some text if the file does not exist
    if (file_in_id == NULL) {
        printf("Not able to open the file.");
    }
    
    line = 0;
    max_x = -10000000;
    max_y = -10000000;
    max_z = -10000000;
    min_x = 10000000;
    min_y = 10000000;
    min_z = 10000000;
    // Read values from the file line by line and store and set min max values for normalization
    while (fscanf(file_in_id, "%f %f %f", &x_file, &y_file, &z_file) == 3) {
        i = line % r;
        j = line / c;

        x[i][j] = x_file;
        if (max_x < x[i][j]) {max_x = x[i][j];}
        if (min_x > x[i][j]) {min_x = x[i][j];}
        y[i][j] = y_file;
        if (max_y < y[i][j]) {max_y = y[i][j];}
        if (min_y > y[i][j]) {min_y = y[i][j];}
        z[i][j] = z_file;
        if (max_z < z[i][j]) {max_z = z[i][j];}
        if (min_z > z[i][j]) {min_z = z[i][j];}
        line++;
    }

    // Close the file
    fclose(file_in_id);
    
    // Normalize x y z grids and init grids
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            // Normalize
            x[i][j] = (x[i][j] - min_x) / (max_x - min_x);
            y[i][j] = (y[i][j] - min_y) / (max_y - min_y);
            z[i][j] = (z[i][j] - min_z) / (max_z - min_z);
            // Init
            h_water[i][j] = rain;
        }
    }
    
    dx = x[0][1] - x[0][0];
    dy = y[1][0] - y[0][0];
    
    iter = 0;
    while (iter < iter_max) {
        //printf("Iteration: %d\n", iter);
        //Dinf_step(r, c, dx, rain, z, h_water, power_coef, erosion_strength, deposition_strength, diffusion_coef);
        D8_step(r, c, dx, rain, z, h_water, power_coef, erosion_strength, deposition_strength, diffusion_coef);
        iter++;  
    } 

    // Unormalize x y z grids
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            x[i][j] = x[i][j]  * (max_x - min_x) + min_x;
            y[i][j] = y[i][j]  * (max_y - min_y) + min_y;
            z[i][j] = z[i][j]  * (max_z - min_z) + min_z;
        }
    }

    write_xyz_to_binary_file(r, c, x, y, z, file_out);

    // End the clock
    clock_t end_time = clock();

    // Calculate the execution time in seconds
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Print the execution time
    printf("Time taken: %.6f seconds\n", time_taken);

    return 0;
}
