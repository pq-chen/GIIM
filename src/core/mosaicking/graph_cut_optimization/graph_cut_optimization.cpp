#include "graph_cut_optimization.h"

#include <cstring>


void Energy::AddDualTerm(int x, int y, int e00, int e01, int e10, int e11) {
  this->add_tweights(x, e11, e00);
  e01 -= e00;
  e10 -= e11;

  if (e01 < 0) {
	  this->add_tweights(x, 0, e01);
	  this->add_tweights(y, 0, -e01);
	  this->add_edge(x, y, 0, e01 + e10);
  } else if (e10 < 0) {
	  this->add_tweights(x, 0, -e10);
	  this->add_tweights(y, 0, e10);
	  this->add_edge(x, y, e01 + e10, 0);
  } else {
	  this->add_edge(x, y, e01, e10);
  }
}

BlockList::~BlockList() {
  while (head_ != 0) {
	auto block(head_);
	head_ = head_->next;
	delete block;
  }
}

void BlockList::AddItem(void* item) {
  if (head_size_ == 4) {
	auto block(new Block);
	block->next = head_;
	head_ = block;
	head_size_ = 0;
  }
  head_->item[head_size_] = item;
  head_size_++;
}

void* BlockList::GetNext() {
  void* item(cursor_->item[cursor_idx_++]);
  if ((cursor_ == head_ && cursor_idx_ >= head_size_) || cursor_idx_ == 4) {
	cursor_ = cursor_->next;
	cursor_idx_ = 0;
  }
  return item;
}

GraphCutOptimization::GraphCutOptimization(int sites_count, int labels_count)
	: labels_count_(labels_count),
	  sites_count_(sites_count),
	  active_sites_index_(new int[sites_count]),
	  labels_(new int[sites_count]),
	  cur_site_costs_(new int[sites_count]) {
  memset(labels_, 0, sites_count * sizeof(int));
  memset(active_sites_index_, -1, sites_count * sizeof(int));
}

GraphCutOptimization::~GraphCutOptimization() {
  delete[] active_sites_index_;
  delete[] labels_;
  delete[] cur_site_costs_;
  if (data_costs_) delete[] data_costs_;
}

long long GraphCutOptimization::Expansion(int iterations_count) {
  UpdateLabelInfo();

  long long new_energy(ComputeEnergy()), old_energy(new_energy + 1);
  for (int i = 1; i <= iterations_count; i++) {
	old_energy = new_energy;
	new_energy = OneExpansionIteration();
	if (new_energy == old_energy) break;
  }

  cur_iteration_ = 0;
  iterations_count_ = 0;
  return new_energy;
}

bool GraphCutOptimization::AlphaExpansion(int alpha_label) {
  FinalizeNeighbors();

  if (iterations_count_ == 0)
	update_label_info_ = true;
  UpdateLabelInfo();

  int* active_sites = new int[sites_count_];
  int active_sites_count(QueryActiveSites(alpha_label, active_sites));
  if (!active_sites_count) {
	delete[] active_sites;
	return false;
  }
  for (int i = 0; i < active_sites_count; i++)
	active_sites_index_[active_sites[i]] = i;

  Energy energy(active_sites_count, total_neighbors_count_);
  energy.AddVars(active_sites_count);
  old_energy_ = 0;
  SetupDataCosts(active_sites_count, alpha_label, &energy, active_sites);
  SetupSmoothCosts(active_sites_count, alpha_label, &energy, active_sites);
  long long new_energy(energy.Minimize());
  if (new_energy < old_energy_)
	ApplyNewLabel(&energy, active_sites, active_sites_count, alpha_label);
  for (int i = 0; i < active_sites_count; i++)
	active_sites_index_[active_sites[i]] = -1;

  delete[] active_sites;
  return new_energy < old_energy_;
}

void GraphCutOptimization::SetDataCost(int site, int lable, int e) {
  if (!data_costs_) {
	int* table = new int[sites_count_ * labels_count_];
	memset(table, 0, sites_count_ * labels_count_ * sizeof(int));
	data_costs_ = table;
	update_label_info_ = true;
  }
  data_costs_[site * labels_count_ + lable] = e;
  if (labels_[site] == lable)
	update_label_info_ = true; // m_labelingDataCosts is dirty
}

void GraphCutOptimization::SetSmoothCost(SmoothCostFunc func, void* data) {
  smooth_cost_func_ = func;
  smooth_data_ = data;
}

int GraphCutOptimization::WhatLabel(int site) {
  return labels_[site];
}

void GraphCutOptimization::SetLabel(int site, int label) {
  labels_[site] = label;
  update_label_info_ = true;
}

long long GraphCutOptimization::ComputeEnergy() {
  return ComputeDataEnergy() + ComputeSmoothEnergy();
}

long long GraphCutOptimization::ComputeDataEnergy() {
  UpdateLabelInfo();
  long long energy = 0;
  for (int i = 0; i < sites_count_; i++)
	energy += cur_site_costs_[i];
  return energy;
}

long long GraphCutOptimization::ComputeSmoothEnergy() {
  FinalizeNeighbors();
  long long e(0);
  for (int i = 0; i < sites_count_; i++) {
	int neighbors_count, *neighbors_idx, *neighbors_weight;
	GetNeighborInfo(i, &neighbors_count, &neighbors_idx, &neighbors_weight);
	for (int j = 0; j < neighbors_count; j++) {
	  int site(neighbors_idx[j]);
	  if (site < i)
		e += neighbors_weight[j] * (smooth_cost_func_(
			i, site, labels_[i], labels_[site], smooth_data_));
	}
  }
  return e;
}

int GraphCutOptimization::QueryActiveSites(
	int alpha_label,
	int* active_sites) {
  int count(0);
  for (int i = 0; i < sites_count_; i++)
	if (labels_[i] != alpha_label)
	  active_sites[count++] = i;
  return count;
}

void GraphCutOptimization::SetupDataCosts(
	int size,
	int alpha_label,
	Energy* energy,
	int* active_sites) {
  for (int i = 0; i < size; i++)
	AddUnitTerm(
		energy, i, data_costs_[active_sites[i] * labels_count_ + alpha_label],
		cur_site_costs_[active_sites[i]]);
}

void GraphCutOptimization::SetupSmoothCosts(
	int size,
	int alpha_label,
	Energy* energy,
	int* active_sites) {
  for (int i = size - 1; i >= 0; i--) {
	int site1(active_sites[i]),
		neighbors_count,
		*neighbors_idx,
		*neighbors_weight;
	GetNeighborInfo(
		site1, &neighbors_count, &neighbors_idx, &neighbors_weight);
	for (int j = 0; j < neighbors_count; j++) {
	  int site2(neighbors_idx[j]);
	  if (active_sites_index_[site2] == -1) {
		AddUnitTerm(
			energy, i,
			smooth_cost_func_(
				site1, site2, alpha_label, labels_[site2], smooth_data_),
			smooth_cost_func_(
				site1, site2, labels_[site1], labels_[site2], smooth_data_),
			neighbors_weight[j]);
	  } else if (site2 < site1) {
		AddDualTerm(
			energy, i, active_sites_index_[site2],
			smooth_cost_func_(
				site1, site2, alpha_label, alpha_label, smooth_data_),
			smooth_cost_func_(
				site1, site2, alpha_label, labels_[site2], smooth_data_),
			smooth_cost_func_(
				site1, site2, labels_[site1], alpha_label, smooth_data_),
			smooth_cost_func_(
				site1, site2, labels_[site1], labels_[site2], smooth_data_),
			neighbors_weight[j]);
	  }
	}
  }
}

void GraphCutOptimization::ApplyNewLabel(
	Energy* energy,
	int* active_sites,
	int size,
	int alpha_label) {
  for (int i = 0; i < size; i++)
	if (energy->GetVar(i) == 0) {
	  int site(active_sites[i]);
	  labels_[site] = alpha_label;
	  cur_site_costs_[site] = data_costs_[site * labels_count_ + alpha_label];
	}
}

void GraphCutOptimization::UpdateLabelDataCosts() {
  for (int i = 0; i < sites_count_; i++)
	cur_site_costs_[i] = data_costs_[i * labels_count_ + labels_[i]];
}

void GraphCutOptimization::UpdateLabelInfo(
	bool updateCounts,
	bool updateActive,
	bool updateCosts) {
  if (!update_label_info_) return;

  update_label_info_ = false;

  if (updateCosts)
	UpdateLabelDataCosts();
}

void GraphCutOptimization::AddUnitTerm(Energy* energy, int i, int e0, int e1) {
  old_energy_ += e1;
  energy->AddUnitTerm(i, e0, e1);
}

void GraphCutOptimization::AddUnitTerm(
	Energy* energy,
	int i,
	int e0,
	int e1,
	int weight) {
  old_energy_ += e1 * weight;
  energy->AddUnitTerm(i, e0 * weight, e1 * weight);
}

void GraphCutOptimization::AddDualTerm(
	Energy* energy,
	int i,
	int j,
	int e00,
	int e01,
	int e10,
	int e11,
	int weight) {
  old_energy_ += e11 * weight;
  energy->AddDualTerm(
	  i, j, e00 * weight, e01 * weight, e10 * weight, e11 * weight);
}

long long GraphCutOptimization::OneExpansionIteration() {
  cur_iteration_ = 0;
  iterations_count_ = labels_count_;
  for (int i = 0; i < labels_count_; i++, cur_iteration_++)
	AlphaExpansion(i);
  return ComputeEnergy();
}

GraphCutOptimizationGeneralGraph::~GraphCutOptimizationGeneralGraph() {
  if (neighbors_)
	delete[] neighbors_;
  if (neighbors_counts_) {
	for (int i = 0; i < sites_count_; i++)
	  if (neighbors_counts_[i] != 0 ) {
		delete[] neighbors_indexes_[i];
		delete[] neighbors_weights_[i];
	  }
	delete[] neighbors_counts_;
	delete[] neighbors_indexes_;
	delete[] neighbors_weights_;
  }
}

void GraphCutOptimizationGeneralGraph::SetNeighbors(
	int site1,
	int site2,
	int weight) {
  if (!neighbors_)
	neighbors_ = new BlockList[sites_count_];
  auto neighbor1(new Neighbor), neighbor2(new Neighbor);
  neighbor1->weight = weight;
  neighbor1->to_node = site2;
  neighbor2->weight = weight;
  neighbor2->to_node = site1;
  neighbors_[site1].AddItem(neighbor1);
  neighbors_[site2].AddItem(neighbor2);
}

void GraphCutOptimizationGeneralGraph::GetNeighborInfo(
	int site,
	int* neighbors_count,
	int** neighbors_idx,
	int** neighbors_weight) {
  if (neighbors_counts_) {
	(*neighbors_count) = neighbors_counts_[site];
	(*neighbors_idx) = neighbors_indexes_[site];
	(*neighbors_weight) = neighbors_weights_[site];
  } else {
	*neighbors_count = 0;
	*neighbors_idx = 0;
	*neighbors_weight = 0;
  }
}

void GraphCutOptimizationGeneralGraph::FinalizeNeighbors() {
  if (!finalize_neighbors_) return;
  finalize_neighbors_ = false;

  int* weights = new int[sites_count_], *idxes = new int[sites_count_];
  neighbors_counts_ = new int[sites_count_];
  neighbors_indexes_ = new int*[sites_count_];
  neighbors_weights_ = new int*[sites_count_];

  for (int site = 0; site < sites_count_; site++) {
	if (neighbors_ && !neighbors_[site].IsEmpty()) {
	  neighbors_[site].SetCursorFront();
	  int count = 0;
	  while (neighbors_[site].HasNext()) {
		auto neighbor(static_cast<Neighbor*>(neighbors_[site].GetNext()));
		idxes[count] = neighbor->to_node;
		weights[count] = neighbor->weight;
		delete neighbor;
		count++;
	  }
	  neighbors_counts_[site] = count;
	  total_neighbors_count_ += count;
	  neighbors_indexes_[site] = new int[count];
	  neighbors_weights_[site] = new int[count];
	  for (int i = 0; i < count; i++) {
		neighbors_indexes_[site][i] = idxes[i];
		neighbors_weights_[site][i] = weights[i];
	  }
	} else {
	  neighbors_counts_[site] = 0;
	}

  }

  delete[] idxes;
  delete[] weights;
  if (neighbors_) {
	delete[] neighbors_;
	neighbors_ = nullptr;
  }
}