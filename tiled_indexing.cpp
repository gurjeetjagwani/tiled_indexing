#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

#define INDEX_4D(N4, N3, N2, N1, I4, I3, I2, I1) \
    (N1 * (N2 * (N3 * I4 + I3) + I2) + I1)
#define INDEX_3D(N3, N2, N1, I3, I2, I1) \
    (N1 * (N2 * I3 + I2) + I1)
#define INDEX_2D(N1, I2, I1) \
    (N1 * I2 + I1)

const int num_pol = 4;

template<typename T>
void prefix_sum(
    const T num_tiles, 
    const T* num_points_in_tiles, 
    T* tile_offsets
)
{
    T sum = 0;
    T i = 0;
    for(i = 0; i < num_tiles; i++)
    {
        T x = num_points_in_tiles[i];
        tile_offsets[i] = sum;
        sum += x;
    }

    tile_offsets[i] = sum;
}

template<typename T, typename F>
static void gen_data(
    T num_baselines, 
    T num_channels, 
    T num_times, 
    std::vector<F> vis,
    std::vector<F> weights, 
    std::vector<F> uvw,
    std::vector<F> freqs
)
{
    for (T i = 0; i < num_times; i++)
    {
        for(T k = 0; k < num_baselines; k++)
        {
            T i_u = INDEX_3D(num_times, num_baselines,3 ,i, k, 0);
            T i_v = INDEX_3D(num_times, num_baselines,3 ,i, k, 1);
                
            for (T l = 0; l < num_channels; l++)
            {
                freqs[l] = 1e9 + l * 1e6;               
                uvw[i_u] = (double)rand() / (double) RAND_MAX;
                uvw[i_v] = (double)rand() / (double) RAND_MAX; 
                
                for (T j = 0; j < num_pol; j++)
                {
                    T i_vis = INDEX_4D(num_times, num_baselines, num_channels, num_pol, i, k, l, j);
                    vis[i_vis] = 1.0;
                    weights[i_vis] = 1.0; 
                }
            }
        }
    }
    return vis;
}
template<typename T, typename F>
static void tile_count_for_indexing( 
    T grid_size, 
    T support,
    T num_channels, 
    T num_baselines,  
    T num_times,
    T top_left_u, 
    T top_left_v,
    T tile_size_u,
    F inv_tile_size_u, 
    F inv_tile_size_v, 
    std::vector<T> num_points_in_tile,
    std::vector<F> uv
)
{
    T grid_centre = grid_size / 2;
    for (T i = 0; i < num_times; i++)
    {
        for(T k = 0; i < num_baselines; k++)
        {
            T i_u = INDEX_3D(num_times, num_baselines,3 ,i, k, 0);
            T i_v = INDEX_3D(num_times, num_baselines,3 ,i, k, 1);
           for(int j=0; j < num_channels; j++)
            { 

                T inv_wavelength = freqs[j] / 299792458.0;  
                T grid_u = round(i_u) * inv_wavelength + grid_centre;
                T grid_v = round(i_v) * inv_wavelength + grid_centre;
            
                if(grid_u + support < grid_size && grid_v + support < grid_size)
                {
                    T min_support_u, min_support_v, max_support_v, max_support_u;
                    T rel_u = grid_u - top_left_u;
                    T rel_v = grid_v - top_left_v;
                    const float u1 = (float) (rel_u - support) * inv_tile_size_u;
                    const float u2 = (float) (rel_v + support + 1) * inv_tile_size_v;
                    const float v1 = (float) (rel_v - support) * inv_tile_size_v;
                    const float v2 = (float) (rel_v + support + 1) * inv_tile_size_v;
                    min_support_u = (int) (floor(u1));
                    min_support_v = (int) (floor(v1));
                    max_support_u = (int) (ceil(u2));
                    max_support_v = (int) (ceil(v2));
                    for (T pv = min_support_v; pv < max_support_v; pv++)
                    {
                        for(T pu = min_support_u; pu < max_support_u; pu++)
                        {
                            num_points_in_tile +=  pu * tile_size_u + pv;  
                        }
                    } 
                }

            }
        }
    }
}

template<typename T, typename F>
static void bucket_sorted_indexing( 
    T grid_size, 
    T support,
    T top_left_u,
    T top_left_v,
    T num_tiles_u,
    const float inv_tile_size_u,
    const float inv_tile_size_v,
    T num_channels, 
    T num_baselines,  
    T num_times,
    std::vector<F> uvw,
    std::vector<F> freqs,
    std::vector<F> vis, 
    std::vector<T> tile_offsets,
    std::vector<T> sorted_vis_index
)
{
    T grid_centre = grid_size / 2;
    for (int i = 0; i < num_times; i++)
    {
        for(int k = 0; i < num_baselines; k++)
        {   
            T i_u = INDEX_3D(num_times, num_baselines,3 ,i, k, 0);
            T i_v = INDEX_3D(num_times, num_baselines,3 ,i, k, 1);
            
            for(int j=0; j < num_channels; j++)
            {
                T inv_wavelength = freqs[j] / 299792458.0;  
                T grid_u = round(i_u) * inv_wavelength + grid_centre; 
                T grid_v = round(i_v) * inv_wavelength + grid_centre;

                for(int l =0; l < num_pol; l++)
                {
                    T i_vis = INDEX_4D(num_times, num_baselines, num_channels, num_pol, i, k, j, l);
                    if(grid_u + support < grid_size && grid_v + support < grid_size)
                    {
                        T min_support_u, min_support_v, max_support_v, max_support_u;
                        T rel_u = grid_u - top_left_u;
                        T rel_v = grid_v - top_left_v;
                        const float u1 = (float) (rel_u - support) * inv_tile_size_u;
                        const float u2 = (float) (rel_v + support + 1) * inv_tile_size_v;
                        const float v1 = (float) (rel_v - support) * inv_tile_size_v;
                        const float v2 = (float) (rel_v + support + 1) * inv_tile_size_v;
                        min_support_u = (T) (floor(u1));
                        min_support_v = (T) (floor(v1));
                        max_support_u = (T) (ceil(u2));
                        max_support_v = (T) (ceil(v2));
                        for (T pv = min_support_v; pv < max_support_v; pv++)
                        {
                            for(T pu = min_support_u; pu < max_support_u; pu++)
                            {
                                T off = tile_offsets[pu + pv * num_tiles_u];
                                tile_offsets[pu + pv * num_tiles_u] += 1;
                                sorted_vis_index[off] = i_vis;
                            }
                        }
                    }

                }
            }
        }
    }
}

int main()
{
    //Define parameters for generating dummy visibilities, weights and uv plane
    int stations = 64;
    int baselines = stations * (stations -1);
    int channels = 30;
    int times = 1;
    std::vector<double> freqs;
    std::vector<double> uvw;
    std::vector<double> vis;
    std::vector<double> weights;
    gen_data<int, double>(baselines,channels,times,vis,weights,uvw,freqs);
    //Define parameters for the grid and tiles
    int grid_size = 1024;
    int grid_centre = grid_centre / 2;
    int tile_size_u = 32;
    int tile_size_v = 16;
    int support = 3;
    int num_tiles_u = grid_size / tile_size_u;
    int num_tiles_v = grid_size / tile_size_v;
    int num_tiles = num_tiles_u * num_tiles_v;
    int c_tile_u = grid_centre / tile_size_u;
    int c_tile_v = grid_centre / tile_size_v;
    // Check if the central tile is in the centre
    int top_left_u = grid_centre - c_tile_u * tile_size_u - tile_size_u/2;
    int top_left_v = grid_centre - c_tile_v * tile_size_v - tile_size_v/2;
    assert(top_left_u <= 0);
    assert(top_left_v <= 0);
    // Calculate the inverse of the tile sizes for fast division later
    const float inv_tile_size_u = 1.0/float(tile_size_u);
    const float inv_tile_size_v = 1.0/float(tile_size_v);
    std::vector<int> num_points_in_tile;
    std::vector<int> sorted_vis_index;
    std::vector<int> tile_offsets;
    tile_count_for_indexing<int, double>(grid_size, support, channels, baselines, times, top_left_u, top_left_v, tile_size_u, inv_tile_size_u, inv_tile_size_v, num_points_in_tile, uvw);
    prefix_sum<int>(num_tiles, num_points_in_tile.data(), tile_offsets.data());
    bucket_sorted_indexing<int, double>(grid_size, support, top_left_u, top_left_v, num_tiles_u, inv_tile_size_u, inv_tile_size_v, channels, baselines, times, uvw, freqs, vis, tile_offsets, sorted_vis_index);
}