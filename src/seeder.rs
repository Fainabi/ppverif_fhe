pub struct IdSeeder {
    global_seed: u128,
}

impl IdSeeder {
    pub fn new(global_seed: u128) -> Self {
        Self {
            global_seed
        }
    }

    pub fn seed(&self, id: u128) -> tfhe::Seed {
        tfhe::Seed(self.global_seed + id)
    }
}
