#include "maxflow/bk_maxflow_v3.04/graph.cpp"
#include "maxflow/bk_maxflow_v3.04/maxflow.cpp"


#define MAX_ENERGYTERM 10000000

typedef Graph<int, int, long long> IntGraph;

class Energy : public IntGraph {
 public:
  Energy(int max_vars_count, int max_edges_count)
	  : IntGraph(max_vars_count, max_edges_count) {}
  Energy(const Energy&) = delete;
  Energy& operator=(const Energy&) = delete;
  ~Energy() = default;

  void AddVars(int count = 1) { IntGraph::add_node(count); }

  void AddUnitTerm(int x, int e0, int e1) { this->add_tweights(x, e1, e0); }

  void AddDualTerm(int x, int y, int e00, int e01, int e10, int e11);

  long long Minimize() { return e_ + IntGraph::maxflow(); }

  int GetVar(int x) { return static_cast<int>(this->what_segment(x)); }

 private:
  long long e_ = 0;
};

class BlockList {
 public:
  BlockList() = default;
  BlockList(const BlockList&) = delete;
  BlockList& operator=(const BlockList&) = delete;
  ~BlockList();

  void AddItem(void* item);

  bool IsEmpty() { return head_ == nullptr; };

  void SetCursorFront() {
	cursor_ = head_;
	cursor_idx_ = 0;
  };

  void* GetNext();

  bool HasNext() { return cursor_ != nullptr; }

 private:
  struct Block {
	void* item[4];
	Block* next;
  };

  Block* head_ = nullptr;
  Block* cursor_;
  char head_size_ = 4;
  char cursor_idx_;
};

class GraphCutOptimization {
 public:
  typedef int(*SmoothCostFunc)(
      int site1, int site2, int lable1, int lable2, void* data);

  GraphCutOptimization(int sites_count, int labels_count);
  GraphCutOptimization(const GraphCutOptimization&) = delete;
  GraphCutOptimization& operator=(const GraphCutOptimization&) = delete;
  virtual ~GraphCutOptimization();

  long long Expansion(int iterations_count);

  bool AlphaExpansion(int alpha_label);

  void SetDataCost(int site, int lable, int e);

  void SetSmoothCost(SmoothCostFunc func, void* data);

  int WhatLabel(int site);

  void SetLabel(int site, int lable);

  long long ComputeEnergy();

  long long ComputeDataEnergy();

  long long ComputeSmoothEnergy();

 protected:
  int QueryActiveSites(int alpha_label, int* active_sites);

  virtual void GetNeighborInfo(
	  int site,
	  int* neighbors_count,
	  int** neighbors_idx,
	  int** neighbors_weight) = 0;
  
  virtual void FinalizeNeighbors() = 0;

  void SetupDataCosts(
	  int size,
	  int alpha_label,
	  Energy* energy,
	  int* active_sites);

  void SetupSmoothCosts(
	  int size,
	  int alpha_label,
	  Energy* energy,
	  int* active_sites);
  
  void ApplyNewLabel(
	  Energy* energy,
	  int* active_sites,
	  int size,
	  int alpha_label);
  
  void UpdateLabelDataCosts();
  
  void UpdateLabelInfo(
	  bool updateCounts = true,
	  bool updateActive = true,
	  bool updateCosts = true);

  void AddUnitTerm(Energy* energy, int i, int e0, int e1);
  
  void AddUnitTerm(Energy* energy, int i, int e0, int e1, int weight);
  
  void AddDualTerm(
	  Energy* energy,
	  int i,
	  int j,
	  int e00,
	  int e01,
	  int e10,
	  int e11,
	  int weight);

  int sites_count_;
  int labels_count_;
  int total_neighbors_count_ = 0;
  int cur_iteration_ = 0;
  int iterations_count_ = 0;
  long long old_energy_;
  int* labels_;
  int* active_sites_index_;
  int* data_costs_ = nullptr;
  int* cur_site_costs_;
  int* neighbors_counts_ = nullptr;

  bool update_label_info_ = true;

  SmoothCostFunc smooth_cost_func_;
  void* smooth_data_;

 private:
  long long OneExpansionIteration();
};

class GraphCutOptimizationGeneralGraph : public GraphCutOptimization {
 public:
  GraphCutOptimizationGeneralGraph(int sites_count, int labels_count)
	  : GraphCutOptimization(sites_count, labels_count) {}
  GraphCutOptimizationGeneralGraph(
	  const GraphCutOptimizationGeneralGraph&) = delete;
  GraphCutOptimizationGeneralGraph& operator=(
	  const GraphCutOptimizationGeneralGraph&) = delete;
  ~GraphCutOptimizationGeneralGraph() override;
           
  void SetNeighbors(int site1, int site2, int weight = 1);

 protected: 
  void GetNeighborInfo(
	  int site,
	  int* neighbors_count,
	  int** neighbors_idx,
	  int** neighbors_weight) override;

  void FinalizeNeighbors() override;

 private:
  struct Neighbor {
	int to_node;
	int weight;
  };
  
  BlockList* neighbors_ = nullptr;
  int** neighbors_indexes_ = nullptr;
  int** neighbors_weights_ = nullptr;
  bool finalize_neighbors_ = true;
};