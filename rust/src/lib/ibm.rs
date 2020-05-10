#![allow(non_snake_case)]

use rand::distributions::{Distribution};
use rand::distributions::uniform::Uniform;
use rand_distr::Exp;
use rand_distr::{Gamma};

use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use std::collections::{BTreeSet, BTreeMap};
use std::cmp::Ordering;
use rand::Rng;
use std::f64::INFINITY;

use indoc::indoc;
use unindent::unindent;

use std::convert::{TryInto, TryFrom};

const INDIVIDUALS_SQL: &str = "INSERT INTO Individuals VALUES (?,?,?,?);";
const INFECTIONS_SQL: &str = "INSERT INTO Infections VALUES (?,?,?);";
const TRANSITIONS_SQL: &str = "INSERT INTO Transitions VALUES (?,?,?,?);";
const RT_INSERT_SQL: &str = indoc!("
    INSERT OR IGNORE INTO RtSufficientStatistics VALUES (?, 0, 0);
");
const RT_INCREMENT_PRIMARY: &str = indoc!("
    UPDATE RtSufficientStatistics SET n_primary = n_primary + 1 WHERE time_discrete = ?;
");
const RT_INCREMENT_SECONDARY: &str = indoc!("
    UPDATE RtSufficientStatistics SET n_secondary = n_secondary + 1 WHERE time_discrete = ?;
");

pub fn to_i64(x: usize) -> i64  {
    x.try_into().unwrap()
}

pub fn weights_to_cdf(w: &Vec<f64>) -> Vec<f64> {
    cumulative_sum(&weights_to_probabilities(w))
}

pub fn weights_to_probabilities(w: &Vec<f64>) -> Vec<f64> {
    let s: f64 = w.iter().sum();
    w.iter().map(|wi| wi / s).collect()
}

pub fn cumulative_sum(v: &Vec<f64>) -> Vec<f64> {
    let mut cs = Vec::with_capacity(v.len());
    for x in v {
        if cs.len() == 0 {
            cs.push(*x);
        }
        else {
            cs.push(*x + cs[cs.len() - 1]);
        }
    }
    cs
}

fn draw_categorical(rng: &mut Xoshiro256PlusPlus, n: usize, cdf: &Vec<f64>) -> usize {
    assert!(cdf.len() == n || cdf.len() == (n - 1));
    
    if n == 1 {
        return 0;
    }
    else {
        let u: f64 = rng.gen();
        for i in 0..(n - 1) {
            if u < cdf[i] {
                return i;
            }
        }
        return n - 1;
    }
}

#[derive(Debug, Clone)]
pub struct VecSet<T> {
    vec: Vec<T>,
    index_map: BTreeMap<T, usize>,
}

impl<T> VecSet<T> where T: Copy + Ord {
    pub fn new() -> Self {
        Self {
            vec: Vec::new(),
            index_map: BTreeMap::new(),
        }
    }
    
    pub fn add(&mut self, item: T) {
        self.index_map.insert(item, self.vec.len());
        self.vec.push(item);
    }
    
    pub fn remove(&mut self, item: T) {
        let index = self.index_map[&item];
    
        // Remove the last item in vec
        let last_item = self.vec.pop().unwrap();
        if index != self.vec.len() {
            // Move it to where the removed item was, and update the index map
            self.vec[index] = last_item;
            self.index_map.insert(last_item, index);
        }
    }
    
    pub fn len(&self) -> usize {
        self.vec.len()
    }
    
    pub fn sample(&self, mut rng: &mut Xoshiro256PlusPlus) -> T {
        self.vec[Uniform::new(0, self.vec.len()).sample(&mut rng)].clone()
    }
}

#[derive(Debug, Copy, Clone)]
struct Event {
    t: f64,
    individual_id: usize,
}

impl Event {
    pub fn new(t: f64, individual_id: usize) -> Self {
        Self { t, individual_id }
    }
}

impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.t == other.t && self.individual_id == other.individual_id
    }
}

impl Eq for Event { }

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.t < other.t {
            Ordering::Less
        }
        else if self.t > other.t {
            Ordering::Greater
        }
        else {
            self.individual_id.cmp(&other.individual_id)
        }
    }
}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone)]
pub struct State {
    pub id: usize,
    pub name: String,
    pub detail: StateDetail,
}

impl State {
    pub fn new_susceptible(id: usize, name: String) -> Self {
        Self { id, name, detail: StateDetail::Susceptible }
    }
    
    pub fn new_final(id: usize, name: String) -> Self {
        Self { id, name, detail: StateDetail::Final }
    }
    
    pub fn new_infected(
        id: usize, name: String,
    ) -> Self {
        Self {
            id,
            name,
            detail: StateDetail::Infected(None),
        }
    }
    
    pub fn is_susceptible(&self) -> bool {
        match self.detail {
            StateDetail::Susceptible => true,
            _ => false,
        }
    }
    
    pub fn is_final(&self) -> bool {
        match self.detail {
            StateDetail::Final => true,
            _ => false,
        }
    }
    
    pub fn is_infected(&self) -> bool {
        match self.detail {
            StateDetail::Infected(_) => true,
            _ => false,
        }
    }
    
    pub fn is_infectious(&self) -> bool {
        if let StateDetail::Infected(Some(infected_state)) = &self.detail {
            infected_state.infectious
        }
        else {
            false
        }
    }
}

#[derive(Debug, Clone)]
pub enum StateDetail {
    Susceptible,
    Final,
    Infected(Option<InfectedState>)
}

#[derive(Debug, Clone)]
pub struct InfectedState {
    pub infectious: bool,
    pub mean_duration: f64,
    pub gamma_shape: f64,
    pub next_state_ids: Vec<usize>,
    pub transition_cdfs: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct Counts {
    _total: usize,
    _total_by_state: Vec<usize>,
    _total_by_ageclass: Vec<usize>,
    counts: Vec<Vec<usize>>,
}

impl Counts {
    pub fn new(n_states: usize, n_ageclasses: usize) -> Self {
        let counts = std::iter::repeat(
            std::iter::repeat(0).take(n_ageclasses).collect()
        ).take(n_states).collect();
        
        Counts {
            _total: 0,
            _total_by_state: std::iter::repeat(0).take(n_states).collect(),
            _total_by_ageclass: std::iter::repeat(0).take(n_ageclasses).collect(),
            counts
        }
    }
    
    pub fn get(&self, state: usize, ageclass: usize) -> usize {
        self.counts[state][ageclass]
    }
    
    pub fn increment(&mut self, state: usize, ageclass: usize, amount: usize) {
        self.counts[state][ageclass] += amount;
        self._total += amount;
        self._total_by_state[state] += amount;
        self._total_by_ageclass[ageclass] += amount;
    }
    
    pub fn decrement(&mut self, state: usize, ageclass: usize, amount: usize) {
        self.counts[state][ageclass] -= amount;
        self._total -= amount;
        self._total_by_state[state] -= amount;
        self._total_by_ageclass[ageclass] -= amount;
    }
    
    pub fn transition(&mut self, from_state: usize, to_state: usize, ageclass: usize) {
        self.counts[from_state][ageclass] -= 1;
        self.counts[to_state][ageclass] += 1;
        self._total_by_state[from_state] -= 1;
        self._total_by_state[to_state] += 1;
    }
    
    pub fn total(&self) -> usize {
        self._total
    }
    
    pub fn total_for_state(&self, state: usize) -> usize {
        self._total_by_state[state]
    }
    
    pub fn total_for_ageclass(&self, ageclass: usize) -> usize {
        self._total_by_ageclass[ageclass]
    }
}

#[derive(Debug, Clone)]
pub struct CIOverN {
    C: Vec<Vec<f64>>,
    I: Vec<f64>,
    N: Vec<f64>,
//    _val: Vec<Vec<f64>>,
//    _row_sums: Vec<f64>,
}

impl CIOverN {
    fn new(C: Vec<Vec<f64>>, I: Vec<usize>, N: Vec<usize>) -> Self {
        let obj = Self {
            C,
            I: I.iter().map(|x| *x as f64).collect(),
            N: N.iter().map(|x| *x as f64).collect(),
        };
        
        obj
    }
    
    fn increment(&mut self, ageclass: usize) {
        self.I[ageclass] += 1.0;
    }
    
    fn decrement(&mut self, ageclass: usize){
        self.I[ageclass] -= 1.0;
    }
    
    fn update_C(&mut self, C: Vec<Vec<f64>>) {
        self.C = C;
    }
    
    fn n_ageclasses(&self) -> usize {
        self.C.len()
    }
    
    fn value(&self, row: usize, col: usize) -> f64 {
        self.C[row][col] * self.I[col] / self.N[col]
    }
    
    fn row(&self, ageclass: usize) -> Vec<f64> {
        (0..self.n_ageclasses()).map(|col| self.value(ageclass, col)).collect()
    }
    
    fn row_sum(&self, ageclass: usize) -> f64 {
        self.row(ageclass).iter().sum()
    }
    
    fn row_cdf(&self, ageclass: usize) -> Vec<f64> {
        let row = self.row(ageclass);
        weights_to_cdf(&row)
    }
}



#[derive(Debug, Copy, Clone)]
struct Individual {
    id: usize,
    ageclass: usize,
    state_id: usize,
    t_infected: Option<f64>,
}

impl Individual {
    fn new(id: usize, ageclass: usize, state_id: usize, t_infected: Option<f64>) -> Self {
        Individual { id, ageclass, state_id, t_infected }
    }
    
    fn update_state(&self, state_id: usize) -> Self {
        Self { state_id, ..*self }
    }
}

pub struct Simulation {
    n_ageclasses: usize,
    states: Vec<State>,
    susceptible_state_id: usize,
    initial_infected_state_id: usize,
    t_change: Vec<f64>,
    beta: Vec<f64>,
    C: Vec<Vec<Vec<f64>>>,
    intervention_index: usize,
    counts: Counts,
    C_I_over_N: CIOverN,
    pub t: f64,
    next_id: usize,
    individuals: BTreeMap<usize, Individual>,
    infectious_individuals: Vec<VecSet<usize>>,
    t_contact: Vec<f64>,
    event_queue: BTreeSet<Event>,
    rng: Xoshiro256PlusPlus,
}

impl Simulation {
    pub fn new(
        n_ageclasses: usize,
        states: Vec<State>,
        susceptible_state_id: usize,
        initial_infected_state_id: usize,
        t_change: Vec<f64>,
        beta: Vec<f64>,
        C: Vec<Vec<Vec<f64>>>,
        initial_counts: Counts,
        rng_seed_opt: Option<u32>,
        db_transaction: &mut rusqlite::Transaction,
        record_all_events: bool,
    ) -> Self {
        let rng_seed = if let Some(rng_seed) = rng_seed_opt {
            rng_seed
        }
        else {
            rand::thread_rng().gen()
        };
    
        db_transaction.execute_batch(&unindent("
            CREATE TABLE Meta (key, value);
            CREATE TABLE Individuals (time REAL, id INTEGER, ageclass INTEGER, initial_state TEXT);
            CREATE TABLE Infections (time REAL, infected_id INTEGER, infectious_id INTEGER);
            CREATE TABLE Transitions (time REAL, id INTEGER, start_state TEXT, end_state TEXT);
            CREATE TABLE Counts (time REAL, state TEXT, ageclass INTEGER, count INTEGER);
            CREATE TABLE RtSufficientStatistics (
                time_discrete INTEGER NOT NULL PRIMARY KEY, n_primary INTEGER, n_secondary INTEGER
            );
        ")).unwrap();
    
        db_transaction.execute(
            "INSERT INTO Meta VALUES ('rng_seed', ?);",
            rusqlite::params![rng_seed]
        ).unwrap();
        
        let n_states = states.len();
        
        let infectious_individuals =
            std::iter::repeat(VecSet::new()).take(n_ageclasses).collect();
        
        let C_I_over_N = CIOverN::new(
            C[0].clone(),
            std::iter::repeat(0).take(n_ageclasses).collect(),
            initial_counts._total_by_ageclass.clone()
        );
        
        let mut sim = Self {
            n_ageclasses,
            states,
            susceptible_state_id,
            initial_infected_state_id,
            t_change,
            beta,
            C,
            intervention_index: 0,
            counts: Counts::new(n_states, n_ageclasses),
            C_I_over_N,
            t: 0.0,
            next_id: 1,
            individuals: BTreeMap::new(),
            infectious_individuals,
            t_contact: std::iter::repeat(INFINITY).take(n_ageclasses).collect(),
            event_queue: BTreeSet::new(),
            rng: Xoshiro256PlusPlus::seed_from_u64(rng_seed as u64),
        };
        
        // Initialize initial infecteds
        sim.initialize_individuals(
            &initial_counts,
            db_transaction,
            record_all_events,
        );
    
        sim.update_contact();
    
        sim
    }
    
    pub fn write_counts(&self, db_transaction: &mut rusqlite::Transaction) {
        for state in &self.states {
            for ageclass in 0..self.n_ageclasses {
                db_transaction.execute(
                    "INSERT INTO Counts VALUES (?, ?, ?, ?);",
                    rusqlite::params![
                        self.t, state.name,
                        to_i64(ageclass + 1),
                        to_i64(self.counts.get(state.id, ageclass))
                    ]
                ).unwrap();
            }
        }
    }
    
    fn initialize_individuals(
        &mut self, initial_counts: &Counts,
        db_transaction: &mut rusqlite::Transaction,
        record_all_events: bool,
    ) {
        let mut insert_individual = db_transaction.prepare(INDIVIDUALS_SQL).unwrap();
        
        for state in self.states.clone() {
            if state.is_susceptible() || state.is_final() {
                for ageclass in 0..self.n_ageclasses {
                    self.counts.increment(state.id, ageclass, initial_counts.get(state.id, ageclass));
                }
            }
            else if state.is_infected() {
                for ageclass in 0..self.n_ageclasses {
                    for _ in 0..initial_counts.get(state.id, ageclass) {
                        self.add_individual(
                            ageclass, &state, true,
                            &mut insert_individual,
                            record_all_events,
                        );
                    }
                }
            }
        }
    }
    
    fn add_individual(
        &mut self, ageclass: usize, state: &State, is_initial: bool,
        insert_individual: &mut rusqlite::Statement,
        record_all_events: bool,
    ) -> usize {
        let id = self.next_id;
        self.next_id += 1;

        let individual = Individual::new(
            id, ageclass, state.id,
            if is_initial { None } else { Some(self.t) }
        );
        self.individuals.insert(id, individual);
        
        if record_all_events {
            insert_individual.execute(
                rusqlite::params![self.t, i64::try_from(id).unwrap(), i64::try_from(ageclass).unwrap(), state.name]
            ).unwrap();
        }

        if state.is_infectious() {
            self.infectious_individuals[ageclass].add(id);
            self.C_I_over_N.increment(ageclass);
        }
        self.insert_transition_event(state, individual.id);
        
        if is_initial {
            self.counts.increment(state.id, ageclass, 1);
        }
        else {
            self.counts.transition(self.susceptible_state_id, state.id, ageclass);
        }

        id
    }
    
    fn insert_transition_event(&mut self, state: &State, individual_id: usize) {
        let t = self.draw_transition_time(state);
        self.event_queue.insert(Event::new(t, individual_id));
    }
    
    fn S(&self, ageclass: usize) -> f64 {
        self.counts.get(self.susceptible_state_id, ageclass) as f64
    }
    
    fn update_contact(&mut self) {
//        println!("update_contact");
        
        for i in 0..self.n_ageclasses {
            let rate = self.beta[self.intervention_index] * self.S(i) * self.C_I_over_N.row_sum(i);
            self.t_contact[i] = self.t + self.draw_exponential(rate);
        }
    }
    
    fn draw_exponential(&mut self, rate: f64) -> f64 {
        assert!(rate >= 0.0);
        if rate == 0.0 {
            INFINITY
        }
        else {
            Exp::new(rate).unwrap().sample(&mut self.rng)
        }
    }
    
    fn draw_gamma(&mut self, shape: f64, scale: f64) -> f64 {
        assert!(shape > 0.0);
        assert!(scale > 0.0);
        Gamma::new(shape, scale).unwrap().sample(&mut self.rng)
    }
    
    fn draw_transition_time(&mut self, state: &State) -> f64 {
        match &state.detail {
            StateDetail::Infected(Some(infected_state)) => {
                self.t + self.draw_gamma(
                    infected_state.gamma_shape,
                    infected_state.mean_duration / infected_state.gamma_shape
                )
            },
            _ => {
                panic!()
            }
        }
    }
    
    fn get_next_contact(&self) -> (f64, Option<usize>) {
        let mut t = INFINITY;
        let mut ageclass_opt = None;
        for i in 0..self.n_ageclasses {
            if self.t_contact[i] < t {
                t = self.t_contact[i];
                ageclass_opt = Some(i);
            }
        }
        (t, ageclass_opt)
    }
    
    pub fn simulate(
        &mut self, t_until: f64, db_transaction: &mut rusqlite::Transaction, record_all_events: bool,
    ) -> bool {
        let mut insert_individual = db_transaction.prepare(INDIVIDUALS_SQL).unwrap();
        let mut insert_infection = db_transaction.prepare(INFECTIONS_SQL).unwrap();
        let mut insert_transition = db_transaction.prepare(TRANSITIONS_SQL).unwrap();
        let mut insert_rt = db_transaction.prepare(RT_INSERT_SQL).unwrap();
        let mut increment_rt_primary = db_transaction.prepare(RT_INCREMENT_PRIMARY).unwrap();
        let mut increment_rt_secondary = db_transaction.prepare(RT_INCREMENT_SECONDARY).unwrap();
        
        let mut done = false;
        while self.t < t_until {
            let mut found_event = false;
            
            let (t_contact, ageclass_opt) = self.get_next_contact();
            let t_transition = self.t_next_transition().unwrap_or(INFINITY);
//            println!("t_contact = {}, t_transition = {}", t_contact, t_transition);
            
            if !t_contact.is_finite() && !t_transition.is_finite() {
                done = true;
                self.t = t_until;
                break;
            }
            
            if t_contact < t_transition {
                if t_contact <= t_until {
                    self.do_contact_event(
                        t_contact, ageclass_opt.unwrap(),
                        &mut insert_individual,
                        &mut insert_infection,
                        &mut insert_transition,
                        &mut insert_rt,
                        &mut increment_rt_primary,
                        &mut increment_rt_secondary,
                        record_all_events,
                    );
                    found_event = true;
                }
            }
            else {
                if t_transition.is_finite() && t_transition <= t_until {
                    let event = self.dequeue_next_transition_event().unwrap();
                    self.do_transition_event(event, &mut insert_transition, record_all_events);
                    found_event = true;
                }
            }
            
            if !found_event {
                self.t = t_until;
            }
            
            // Identify intervention changepoint
            if self.intervention_index < self.C.len() - 1 {
                if self.t > self.t_change[self.intervention_index] {
                    self.intervention_index += 1;
                    self.C_I_over_N.update_C(self.C[self.intervention_index].clone());
                    self.update_contact();
                    
                    println!("Updated intervention to {} at t = {}", self.intervention_index, self.t);
                }
            }
        }
    
        done
    }
    
    pub fn do_contact_event(
        &mut self, t: f64, ageclass: usize,
        insert_individual: &mut rusqlite::Statement,
        insert_infection: &mut rusqlite::Statement,
        insert_transition: &mut rusqlite::Statement,
        insert_rt: &mut rusqlite::Statement,
        increment_rt_primary: &mut rusqlite::Statement,
        increment_rt_secondary: &mut rusqlite::Statement,
        record_all_events: bool,
    ) {
//        println!("do_contact_event()");
        self.t = t;
        
        // Draw ageclass of infecting individual proportional to C_I_over_N
        let cdf = self.C_I_over_N.row_cdf(ageclass);
//        println!("cdf: {:?}", cdf);
        let infecting_ageclass = draw_categorical(&mut self.rng, self.n_ageclasses, &cdf);
//        println!("infecting_ageclass = {}", infecting_ageclass);
        
        // Randomly choose an infectious individual from the infecting ageclass
        let infectious_id = self.infectious_individuals[infecting_ageclass].sample(&mut self.rng);
        let infectious_individual = self.individuals[&infectious_id];
        
        // Create a new infected individual
        let state = self.states[self.initial_infected_state_id].clone();
        let infected_id = self.add_individual(
            ageclass, &state, false,
            insert_individual,
            record_all_events
        );
        
        // Insert individual events
        if record_all_events {
            insert_infection.execute(
                rusqlite::params![self.t, i64::try_from(infected_id).unwrap(), i64::try_from(infectious_id).unwrap()]
            ).unwrap();
            insert_transition.execute(
                rusqlite::params![
                    self.t, i64::try_from(infected_id).unwrap(),
                    self.states[self.susceptible_state_id].name,
                    state.name
                ]
            ).unwrap();
        };
        
        // Update count of people infected during this discrete timestep
        // (denominator of Rt)
        let t_discrete_present = self.t.ceil() as i64;
        insert_rt.execute(rusqlite::params![t_discrete_present]).unwrap();
        increment_rt_primary.execute(rusqlite::params![t_discrete_present]).unwrap();
        
        // Update count of number of infections caused by people infected at a previous timestep
        // (numerator of Rt)
        if let Some(t_past) = infectious_individual.t_infected {
            let t_past_discrete = t_past.ceil() as i64;
            increment_rt_secondary.execute(rusqlite::params![t_past_discrete]).unwrap();
        }
        
        // Update contact times
        self.update_contact()
    }
    
    fn do_transition_event(
        &mut self, event: Event,
        insert_transition: &mut rusqlite::Statement,
        record_all_events: bool,
    ) {
//        println!("do_transition_event()");
        
        self.t = event.t;
        
        let id = event.individual_id;
        let individual = self.individuals[&id].clone();
        let ageclass = individual.ageclass;
        let last_state = self.states[individual.state_id].clone();
        let last_infected_state = match last_state.detail {
            StateDetail::Infected(Some(ref infected_state)) => {
                infected_state
            },
            _ => {
                panic!()
            }
        };
        
        let next_state_id = last_infected_state.next_state_ids[
            draw_categorical(
                &mut self.rng,
                last_infected_state.next_state_ids.len(),
                &last_infected_state.transition_cdfs[individual.ageclass]
            )
        ];
//        println!("last state: {}; next state: {}", self.states[last_state.id].name, self.states[next_state_id].name);
    
        let next_state = self.states[next_state_id].clone();
        assert!(next_state.is_infected() || next_state.is_final());
        
        // Update ageclass-specific state counts
//        println!("Transitioning {} to {}", last_state.id, next_state.id);
        self.counts.transition(last_state.id, next_state.id, ageclass);
        
        // Update infectious counts
        match (last_state.is_infectious(), next_state.is_infectious()) {
            (false, true) => {
                self.infectious_individuals[ageclass].add(id);
                self.C_I_over_N.increment(ageclass);
            },
            (true, false) => {
                self.infectious_individuals[ageclass].remove(id);
                self.C_I_over_N.decrement(ageclass);
            },
            _ => {},
        }
        
        if next_state.is_final() {
            // Remove individual from memory if they're moving to a final state
            self.individuals.remove(&id);
        }
        else {
            // Otherwise insert an updated individual and queue the next transition
            self.individuals.insert(id, individual.update_state(next_state.id));
            self.insert_transition_event(&next_state, individual.id);
        }
        
        if record_all_events {
            insert_transition.execute(
                rusqlite::params![
                    self.t,
                    i64::try_from(id).unwrap(),
                    last_state.name,
                    next_state.name
                ]
            ).unwrap();
        };
        
        // Update contact times
        self.update_contact()
    }
    
    fn t_next_transition(&self) -> Option<f64> {
        self.event_queue.iter().next().map(|event| event.t)
    }
    
    fn dequeue_next_transition_event(&mut self) -> Option<Event> {
        let event_opt = self.event_queue.iter().next().map(|event| *event);
        if let Some(event) = event_opt {
            self.event_queue.remove(&event);
        }
        event_opt
    }
}

#[cfg(test)]
mod tests {
    use crate::ibm::{weights_to_cdf, draw_categorical};
    use rand_xoshiro::Xoshiro256PlusPlus;
    use rand_xoshiro::rand_core::SeedableRng;
    
    #[test]
    fn test_draw_categorical() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        
        let weights = vec![4.0, 3.0, 2.0, 1.0];
        println!("weights: {:?}", weights);
        let cdf = weights_to_cdf(&weights);
        println!("cdf: {:?}", cdf);
        
        let mut counts: Vec<usize> = vec![0, 0, 0, 0];
        for _i in 0..100000 {
            counts[draw_categorical(&mut rng, 4, &cdf)] += 1;
        }
        println!("counts: {:?}", counts);
    }
}
