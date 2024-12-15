
#include <stdlib.h>
#include <stdio.h>
#include <stdalign.h> // Required for alignof in C11
#include <time.h>
#include <math.h>

size_t round_up_to( size_t n, size_t multiple ) {
  size_t const remainder = n % multiple;
  return remainder == 0 ? n : n + multiple - remainder;
}

void** matrix2d_new( size_t esize, size_t ealign, size_t idim, size_t jdim ) {
  // ensure &elements[0][0] is suitably aligned
  size_t const ptrs_size = round_up_to( sizeof(void*) * idim, ealign );
  size_t const row_size = esize * jdim;
  // allocate the row pointers followed by the elements
  void **const rows = malloc( ptrs_size + idim * row_size );
  char *const elements = (char*)rows + ptrs_size;
  for ( size_t i = 0; i < idim; ++i )
    rows[i] = &elements[ i * row_size ];
  return rows;
}

int is_within_bounds(int i, int j, int rows, int cols) {
    return (i >= 0 && i< rows && j >= 0 && j < cols);
}


void D8_update(int r, int c, float dx, float dy, float rain, float z[r][c], int d8_direction_map[r][c], float d8_slope_map[r][c], float h_water[r][c]) {
  int i,j;
  int neighbour_i, neighbour_j, neighbour_id;
  int max_id, max_i, max_j;
  int max_neighbour_id;
  float max_slope, neighbour_slope;
  float d8_distances[8] = {dy, sqrt(dx*dx + dy*dy), dx, sqrt(dx*dx + dy*dy), dy, sqrt(dx*dx + dy*dy), dx, sqrt(dx*dx + dy*dy)};
  int d8_directions[8][2] = {{-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
  float h_water_tampon[r][c];

  // Set grids
  for (i = 0; i < r; i++) {
      for (j = 0; j < c; j++) {
          h_water_tampon[i][j] = h_water[i][j] + rain;
          h_water[i][j] = 0;
          d8_direction_map[i][j] = -1;
          d8_slope_map[i][j] = -1;
      }
  }

  // Update the grid
  for (i = 0; i < r; i++) {
      for (j = 0; j < c; j++) {
          if (h_water_tampon[i][j] > 0) {
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
                          max_id = neighbour_id;
                          max_i = neighbour_i;
                          max_j = neighbour_j;
                      }
                  }
              }
              if (max_slope > 0) {
                  d8_direction_map[i][j] = max_id;
                  d8_slope_map[i][j] = max_slope;
                  h_water[max_i][max_j] += h_water_tampon[i][j];
              }
          }
      }
  }
}

void update_topo(int r, int c, float dx, float dy, float power_coef, float erosion_strength, float deposition_strength, float diffusion_coef, float z[r][c], int d8_direction_map[r][c], float d8_slope_map[r][c], float h_water[r][c]) {
  int i,j;
  for (i = 1; i < r-1; i++) {
      for (j = 1; j < c-1; j++) {
          z[i][j] += diffusion_coef*(z[i][j+1] - 2*z[i][j] + z[i][j-1])/(pow(dx,2)) + diffusion_coef*(z[i+1][j] - 2*z[i][j] + z[i-1][j])/pow(dy,2);
          if (d8_direction_map[i][j] <= 0 && h_water[i][j] > 0) {
              z[i][j] += deposition_strength ;
          }
          else if (d8_direction_map[i][j] > 0 && h_water[i][j] > 0) {
              z[i][j] -=  pow(d8_slope_map[i][j],power_coef)* h_water[i][j]*erosion_strength;
          }
      }
  }
} 
 
    
int main() {
    
    // InOut
    char file_in[] = "file.xyz", file_out[] = "file_out.xyz"; // Change this to your file name
    
    // Simulation parameters
    int  iter, iter_max; // Control iterations
    float deposition_scale, erosion_scale, rain, diffusion_coef, power_coef, erosion_strength, deposition_strength; // Deposition erosion and rain
    
    // Topo grid parameters
    int r = 300, c = 300; // Number of rows and columns
    int i, j; // Loop counters
    int line; // Read file line number
    float x_file, y_file, z_file; // Store fscanf values
    float max_x, max_y, max_z, min_x, min_y, min_z; // Store min/max values
    float dx, dy; // Store dx and dy
    float x[r][c], y[r][c], z[r][c]; // x, y and z grids
    float h_water[r][c]; // water height grid

    // D8 algorithm parameters
    int d8_direction_map[r][c]; //D8 direction map
    float d8_slope_map[r][c]; // D8 slope map
    
    // Initialize parameters 
    iter_max = 500;
    rain = 0.00001;
    diffusion_coef = 0.00000001;
    power_coef = 0.5;
    erosion_strength = 0.1;
    deposition_strength = 0.001;

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
    
    // Normalize x y z grids and init grids
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            // Normalize
            x[i][j] = (x[i][j] - min_x) / (max_x - min_x);
            y[i][j] = (y[i][j] - min_y) / (max_y - min_y);
            z[i][j] = (z[i][j] - min_z) / (max_z - min_z);
            // Init
            h_water[i][j] = rain;
            d8_direction_map[i][j] = -1;
            d8_slope_map[i][j] = -1;
        }
    }
    
    dx = x[0][1] - x[0][0];
    dy = y[1][0] - y[0][0];

    // Close the file
    fclose(file_in_id);
    
    iter = 0;
    while (iter < iter_max) {
        //printf("Iteration: %d\n", iter);
        D8_update(r, c, dx, dy, rain, z, d8_direction_map, d8_slope_map, h_water);
        update_topo(r, c, dx, dy, power_coef, erosion_strength, deposition_strength, diffusion_coef, z, d8_direction_map, d8_slope_map, h_water);
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

    FILE *file_out_id = fopen(file_out, "w");

    if (file_out_id == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    // Write the X, Y, Z grid values to the file
    for (j = 0; j < c; j++) {
        for (i = 0; i < r; i++) {
            // Write X, Y, Z coordinates to the file
            fprintf(file_out_id, "%f %f %f\n", x[i][j], y[i][j], z[i][j]);
        }
    }

    // Close the file
    fclose(file_out_id);

    printf("Data written to output.xyz\n");

    // End the clock
    clock_t end_time = clock();

    // Calculate the execution time in seconds
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Print the execution time
    printf("Time taken: %.6f seconds\n", time_taken);

    return 0;
}
