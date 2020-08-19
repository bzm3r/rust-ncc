#[derive(Clone)]
pub struct Mat<T: Copy + Default> {
    pub n_r: usize,
    pub n_c: usize,
    dat: Vec<T>,
}

impl<T: Copy + Default> Mat<T> {
    pub fn new(n_r: usize, n_c: usize, default: T) -> Mat<T> {
        Mat {
            n_r,
            n_c,
            dat: vec![default; n_r * n_c],
        }
    }

    fn calc_ix(&self, i: usize, j: usize) -> usize {
        i * self.n_c + j
    }

    pub fn set(&mut self, i: usize, j: usize, x: T) {
        let ix = self.calc_ix(i, j);
        self.dat[ix] = x;
    }

    pub fn get(&self, i: usize, j: usize) -> T {
        let ix = self.calc_ix(i, j);
        self.dat[ix]
    }
}

#[derive(Clone)]
pub struct SymMat<T: Copy + Default> {
    pub n: usize,
    dat: Vec<T>,
}

impl<T: Copy + Default> SymMat<T> {
    pub fn new(n: usize, default: T) -> SymMat<T> {
        SymMat {
            n,
            dat: vec![default; n + (n * (n - 1) / 2)],
        }
    }

    fn calc_ix(&self, i: usize, j: usize) -> usize {
        let (i, j) = if i < j { (i, j) } else { (j, i) };
        let k = self.n - i;
        let a = k + (k * (k - 1) / 2);
        let b = j - i;
        self.dat.len() - a + b
    }

    pub fn set(&mut self, i: usize, j: usize, x: T) {
        let ix = self.calc_ix(i, j);
        self.dat[ix] = x;
    }

    pub fn get(&self, i: usize, j: usize) -> T {
        let ix = self.calc_ix(i, j);
        self.dat[ix]
    }
}
