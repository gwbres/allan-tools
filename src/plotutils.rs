use gnuplot::{Figure, Caption};
use gnuplot::{Color, PointSymbol, LineStyle, DashType};
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
    fg.save_to_png(filepath, 640, 480).unwrap();
}

/// plots given data into file
pub fn plot2d (datasets: Vec<(&Vec<f64>,&Vec<f64>,&Vec<f64>)>,
        plot_title: &str, captions: Vec<&str>, filepath: &str) 
{
    let mut fg = Figure::new();
    fg.set_title(plot_title);
    for i in 0..datasets.len() {
        fg.axes2d()
            .set_x_grid(true)
            .set_x_log(Some(10.0_f64))
            .set_y_grid(true)
            .set_y_log(Some(10.0_f64))
            .y_error_lines(datasets[i].0, datasets[i].1, datasets[i].2,
            &[
                Caption(captions[i]),
                Color("blue"), 
                PointSymbol('x'),
            ]);
    }
    fg.save_to_png(filepath, 640, 480).unwrap();
}

pub fn plot3corner (
    taus: &Vec<f64>,
        adev_ab: (&Vec<f64>,&Vec<f64>),
        data_a: (&Vec<f64>,&Vec<f64>),
            adev_bc: (&Vec<f64>,&Vec<f64>),
            data_b: (&Vec<f64>,&Vec<f64>),
                adev_ca: (&Vec<f64>,&Vec<f64>),
                data_c: (&Vec<f64>,&Vec<f64>))
{
    let mut fg = Figure::new();
    fg.set_title("3 corner");
    fg.axes2d()
        .set_x_grid(true)
        .set_x_log(Some(10.0_f64))
        .set_y_grid(true)
        .set_y_log(Some(10.0_f64))
        .y_error_lines(taus, adev_ab.0, adev_ab.1,
        &[
            Caption("adev(a/b)"),
            Color("blue"), 
            PointSymbol('x'),
        ])
        .y_error_lines(taus, data_a.0, data_a.1,
        &[
            Caption("dev(a)"),
            Color("blue"), 
            PointSymbol('x'),
            LineStyle(DashType::Dash),
        ])
        .y_error_lines(taus, adev_bc.0, adev_bc.1,
        &[
            Caption("adev(b/c)"),
            Color("green"), 
            PointSymbol('x'),
        ])
        .y_error_lines(taus, data_b.0, data_b.1,
        &[
            Caption("dev(b)"),
            Color("green"), 
            PointSymbol('x'),
            LineStyle(DashType::Dash),
        ])
        .y_error_lines(taus, adev_ca.0, adev_ca.1,
        &[
            Caption("adev(c/a)"),
            Color("orange"), 
            PointSymbol('x'),
        ])
        .y_error_lines(taus, data_c.0, data_c.1,
        &[
            Caption("dev(c)"),
            Color("orange"), 
            PointSymbol('x'),
            LineStyle(DashType::Dash),
        ]);
    fg.save_to_png("tests/3corner.png", 640, 480).unwrap();
}
