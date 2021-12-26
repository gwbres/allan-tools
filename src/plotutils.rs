use gnuplot::{Figure, Caption};
use gnuplot::{Color, PointSymbol};
use gnuplot::{PointSize, AxesCommon};

/// plots given data into file
pub fn plot1d (data: Vec<f64>, plot_title: &str, caption: &str, filepath: &str) {
    let mut fg = Figure::new();
    fg.set_title(plot_title);
    let x: Vec<u32> = (0..data.len() as u32).collect();
    fg.axes2d()
        .set_x_grid(true)
        .set_y_grid(true)
        .lines(&x, &data, &[
            Caption(caption),
            Color("blue"), 
            PointSize(10.0),
            PointSymbol('x'),
        ]);
    fg.save_to_png(filepath, 640, 480);
}