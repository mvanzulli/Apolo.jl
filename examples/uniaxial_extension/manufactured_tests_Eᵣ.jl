
using Test
using Apolo
"This functions computes by hand the functional value checking:
- Optical flow hypothesis
- Numerical vs analytical displacements
- Numerical vs analytical reference intensities
- Numerical vs analytical deformed intensities
- Numerical vs analytical functional OpticalFlow value for `E = Eᵣ`
"
function manufactured_tests_Eᵣ(
    intensity_func::Function,
    in_img_func::Function,
    uₗ::Function,
    img_data::ImageData,
    gold_solution::AbstractForwardProblemSolution,
    analytic_f::Function,
)
    # --------------------------
    # Tests problem
    # --------------------------
    # extract reference and deformed Images
    img_ref = reference_img(img_data)
    imgs_def = deformed_imgs(img_data)

    # check optical flow hypothesis
    # --------------------------
    p = (rand(start_roi[1]:Lⱼₛ/20:finish_roi[1]), rand(start_roi[2]:Lⱼₛ/20:finish_roi[2]))
    u_p = (uₗ(p[1], 1), 0.0)
    def_p = p .+ u_p
    @test intensity_func(p..., 0) ≈ intensity_func(def_p..., 1) atol = 1e-6
    p = Tuple(rand.(coords))
    @test img_ref([p]) ≈ [intensity_func(p..., 0)] rtol = 1e-3
    @test imgs_def[1]([p]) ≈ [intensity_func(p..., 1)] rtol = 1e-3

    # get roi coordinates, intensity and displacements
    # ------------------------------------------------
    # cooridnates
    nodes_roi = roi_nodes(img_data)
    roi_vec_coords = roi_nodes_coords(img_data)
    roi_vec_coords_x = getindex.(roi_vec_coords, 1)

    # get displacements of roi_coordinates
    disp_roi = gold_solution(roi_vec_coords)
    disp_roi_numeric_x = getindex.(disp_roi, 1)
    # compare them with analytical solution
    disp_roi_analytic_x = uₗ(roi_vec_coords_x, 1)
    @test disp_roi_analytic_x ≈ disp_roi_numeric_x rtol = 1e-5

    # test reference intensity values
    # --------------------------
    # replace values outside the roi
    int_ref_roi_analytic = [intensity_func(r, 0.5, 0) for r in roi_vec_coords_x]
    # intensity array numeric
    int_ref_roi_numeric = img_ref(roi_vec_coords)
    @test int_ref_roi_analytic ≈ int_ref_roi_numeric rtol = 1e-4

    # test deformed coordinates
    #-------------------------
    # compute deformed coordinates as a vector of tuples
    def_roi_numeric = similar(roi_vec_coords)
    for i in 1:length(disp_roi)
        def_roi_numeric[i] = Tuple(disp_roi[i]) .+ roi_vec_coords[i]
    end

    def_roi_numeric_x = getindex.(def_roi_numeric, 1)
    def_roi_analytic_x = roi_vec_coords_x + disp_roi_analytic_x
    @test def_roi_numeric_x ≈ def_roi_analytic_x atol = 1e-5

    # test deformed intensity values
    # --------------------------
    int_def_roi_analytic = [in_img_func([x, 0.5]) ? intensity_func(x, 0.5, 1) : 0.0 for x in def_roi_analytic_x]
    img_def = imgs_def[1]
    int_def_roi_numeric = img_def(def_roi_numeric)
    @test int_def_roi_numeric ≈ int_def_roi_analytic rtol = 1e-1

    # indices of roi points which ⊂ img and not
    #--------------------------------------------
    not_in_img = findall(x -> !in_img_func([x, 0.5]), def_roi_analytic_x)
    in_img = findall(x -> in_img_func([x, 0.5]), def_roi_analytic_x)

    # test error from to the points which are not inside the image
    # functional of points that not remmain inside the image
    #--------------------------------------------
    not_in_img_error_analytic = sum(int_ref_roi_analytic[not_in_img] .^ 2)
    not_in_img_error_numeric = sum(int_ref_roi_numeric[not_in_img] .^ 2)
    @test not_in_img_error_analytic ≈ not_in_img_error_numeric rtol = 1e-4

    # functional of points that remmain inside the image
    #--------------------------------------------
    msf_analytic_inimg = sum((int_def_roi_analytic[in_img] - int_ref_roi_analytic[in_img]) .^ 2)
    msf_numeric_inimg = sum((int_def_roi_numeric[in_img] - int_ref_roi_numeric[in_img]) .^ 2)

    # functional values
    #--------------------------------------------
    msf_numeric = sum((int_def_roi_numeric - int_ref_roi_numeric) .^ 2)
    msf_analytic = sum((int_def_roi_analytic - int_ref_roi_analytic) .^ 2)

    # tests
    #-------
    # analytic expression of the functional extracted form (https://github.com/jorgepz/Materialis.jl/blob/main/examples/extension/extension.jl)
    tf = 1.0

    @test msf_analytic_inimg ≈ analytic_f(Eᵣ, tf) atol = 1e-8
    @test msf_numeric_inimg ≈ analytic_f(Eᵣ, tf) atol = 1e-8
    @test msf_numeric ≈ msf_analytic atol = 1e-2


    return msf_numeric
end
