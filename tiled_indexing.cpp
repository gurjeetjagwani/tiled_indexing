#include <iostream>
#include <cmath>
#include <vector>

#define INDEX_4D(N4, N3, N2, N1, I4, I3, I2, I1) \
    (N1 * (N2 * (N3 * I4 + I3) + I2) + I1)
#define INDEX_3D(N3, N2, N1, I3, I2, I1) \
    (N1 * (N2 * I3 + I2) + I1)
#define INDEX_2D(N1, I2, I1) \
    (N1 * I2 + I1)

const int num_pol = 4;
const int max_abs_uv = 1e9;

template<typename T, typename F>
static void gen_data(
    T num_baselines, 
    T num_channels, 
    T num_times, 
    std::vector<F> vis,
    std::vector<F> weights, 
    std::vector<F> uv
)
{
    for (T i = 0; i < num_times; i++)
    {
        for(T k = 0; k < num_baselines; k++)
        {
            for (T l = 0; l < num_channels; l++)
            {
                for (T j = 0; j < num_pol; j++)
                {
                    T i_vis = INDEX_4D(num_times, num_baselines, num_channels, num_pol, i, k, l, j)
                    vis[i_vis] = 1.0;
                    weights[i_vis] = 1.0; 
                    T i_uv_u = INDEX_2D(2,i_vis,0);
                    T i_uv_v = INDEX_2D(2,i_vis,1);
                    uv[i_uv_u] = (double)rand() / (double) RAND_MAX;
                    uv[i_uv_v] = (double)rand() / (double) RAND_MAX;
                }
            }
        }
    }
    return vis;
}

template<typename T, typename F>
static void bucket_sorted_indexing( 
    T grid_size, 
    T support,
    F cell_size_rad
)
{
    T grid_centre = grid_size / 2;
    T tile_size_u = 32; 
    T tile_size_v = 16;
    F grid_scale = grid_size * cell_size_rad;

}

int main()
{
    //Define parameters for generating dummy visibilities, weights and uv plane
    int stations = 64;
    int baselines = stations * (stations -1);
    int channels = 128;
    int times = 1;
    std::vector<double> uv;
    std::vector<double> vis;
    std::vector<double> weights;
    gen_data<int, double>(baselines,channels,times,vis,weights,uv);
    //Define parameters for the grid
    int grid_size = 1024;
}