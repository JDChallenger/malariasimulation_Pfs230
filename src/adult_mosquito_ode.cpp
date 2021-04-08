/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "adult_mosquito_ode.h"

AdultMosquitoModel::AdultMosquitoModel(
    MosquitoModel growth_model,
    double mu,
    double tau) : growth_model(growth_model), mu(mu), tau(tau)
{}

integration_function_t create_ode(AdultMosquitoModel& model) {
    //create original ode
    auto growth_ode = create_ode(model.growth_model);
    return [&model, growth_ode](const state_t& x, state_t& dxdt, double t) {
        //set submodel total_M
        Rcpp::Rcerr << "current_total_M " << model.growth_model.total_M << std::endl;
        Rcpp::Rcerr << "state size " << x.size() << std::endl;
        for (auto i : x) {
            Rcpp::Rcerr << "in state " << i << std::endl;
        }
        Rcpp::Rcerr << "indices " << get_idx(AdultODEState::S) << std::endl;
        Rcpp::Rcerr << "indices " << get_idx(AdultODEState::E) << std::endl;
        Rcpp::Rcerr << "indices " << get_idx(AdultODEState::I) << std::endl;
        Rcpp::Rcerr << "values " << x[get_idx(AdultODEState::S)] << std::endl;
        Rcpp::Rcerr << "values " << x[get_idx(AdultODEState::E)] << std::endl;
        Rcpp::Rcerr << "values " << x[get_idx(AdultODEState::I)] << std::endl;
        model.growth_model.total_M =
            x[get_idx(AdultODEState::S)] +
            x[get_idx(AdultODEState::E)] +
            x[get_idx(AdultODEState::I)];

        //run the growth ode
        growth_ode(x, dxdt, t);

        //run the adult ode
        auto incubation_survival = exp(-model.mu * model.tau);
        dxdt[get_idx(AdultODEState::S)] =
            .5 * x[get_idx(ODEState::P)] / model.growth_model.dp //growth to adult female
            - x[get_idx(AdultODEState::S)] * model.foim //infections
            - x[get_idx(AdultODEState::S)] * model.mu; //deaths   
        dxdt[get_idx(AdultODEState::E)] =
            x[get_idx(AdultODEState::S)] * model.foim  //infections
            - model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[get_idx(AdultODEState::E)] * model.mu; // deaths
        dxdt[get_idx(AdultODEState::I)] = model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[get_idx(AdultODEState::I)] * model.mu; // deaths
    };
}

//[[Rcpp::export]]
Rcpp::XPtr<AdultMosquitoModel> create_adult_mosquito_model(
    Rcpp::XPtr<MosquitoModel> growth_model,
    double mu,
    double tau
    ) {
    auto model = new AdultMosquitoModel(*growth_model, mu, tau);
    return Rcpp::XPtr<AdultMosquitoModel>(model, true);
}

//[[Rcpp::export]]
void adult_mosquito_model_update(
    Rcpp::XPtr<AdultMosquitoModel> model,
    double mu,
    double foim,
    double susceptible
    ) {
    model->mu = mu;
    model->foim = foim;
    model->lagged_incubating.pop();
    model->lagged_incubating.push(susceptible * foim);
}

//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_adult_solver(
    Rcpp::XPtr<AdultMosquitoModel> model,
    std::vector<double> init) {
    return Rcpp::XPtr<Solver>(
        new Solver(init, create_ode(*model)),
        true
    );
}
