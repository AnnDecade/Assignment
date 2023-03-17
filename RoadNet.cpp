#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <queue>
#include <utility>
#include <string>
using namespace std;
#define REP(i, n) for(int i=0; i<n; ++i)
#define FOR(i, a, b) for(int i=a; i<b; ++i)
const int INF = 1e9;
const int NONE = -1;

// coordination
struct Point {
	double x, y;

	bool operator<(const Point& e) const {
		if(x == e.x) return y < e.y;
		return x < e.x;
	}
	bool operator>(const Point& e) const {
		if(x == e.x) return y > e.y;
		return x > e.x;
	}
	bool operator==(const Point& e) const {
		return x==e.x && y==e.y;
	}
};

// line segment
struct Segment {
	Point left, right;

	bool operator==(const Segment& e) const {
		return left==e.left && right == e.right;
	}
};

// intersection info
// if flg is false, this represents T-junction
// if flg is false, p.x stores point number
struct Intersection {
	Point p;
	mutable set<int> s; // which lines intersect
	mutable bool flg = true; // whether intersection, or not
	
	bool operator<(const Intersection& e) const {
		return p < e.p;
	}
	bool operator==(const Intersection& e) const {
		return p == e.p;
	}
};

// stored path and its distance
struct PathInfo {
	double cost;
	vector<int> path;

	bool operator<(const PathInfo& e) const {
		if(cost == e.cost) {
			int size = min(path.size(), e.path.size());
			REP(i, size) {
				if(path[i] == e.path[i]) continue;
				return path[i] < e.path[i];
			}
			return path.size() < e.path.size();
		}
		return cost < e.cost;
	}
	bool operator==(const PathInfo& e) const {
		return cost == e.cost && path == e.path;
	}
};

template<typename T>
size_t getIdx(vector<T> v, T n) {
	auto itr = find(v.begin(), v.end(), n);

	return distance(v.begin(), itr);
}

// get Point number for array(if intersection, head word is 'C')
int Point2Num(string str, int N, int c) {
	if(str[0] != 'C') {
		if(stoi(str) <= N) return (stoi(str)-1);
		else return (stoi(str)+c-1);
	}
	else {
		str.erase(0, 1);
		return (stoi(str)+N-1);
	}
}

// get Point name for output
string Num2Point(int num, int N, int c, int P) {
	if(num < N) return to_string(num+1);
	else if(num < N+c) return 'C'+to_string(num-N+1);
	else return to_string(num-c+1);
}

// check head words whether corresponding to key
bool checkHWords(vector<int> vec, vector<int> key) {
       for(int i=0; i<key.size(); ++i)  {
	       if(vec[i] != key[i]) return false;
       }
       return true;
}	       

// This function calculates destance between 2 points
double calDistance(Point p, Point q) {
	double dx = abs(p.x-q.x);
	double dy = abs(p.y-q.y);

	return sqrt(dx*dx + dy*dy);
}

int orientation(Point p, Point q, Point r) {
	int val = (q.y-p.y)*(r.x-q.x)-(q.x-p.x)*(r.y-q.y);
	
	if(val == 0) return 0;    // collinear

	return (val>0) ? 1 : 2;   // CW or CCW
}

// whether intersection or not
//  1 : detect intersection
//  0 : doesn't detect intersection
// negative number : T-junction
int doIntersect(Segment s1, Segment s2) {
	Point p1 = s1.left, q1 = s2.left, p2 = s1.right, q2 = s2.right;

	// definition
	int o1 = orientation(p1, p2, q1);
	int o2 = orientation(p1, p2, q2);
	int o3 = orientation(q1, q2, p1);
	int o4 = orientation(q1, q2, p2);

	if(o1!=o2 && o3!=o4) {
		// doesn't include collinear in any 3 points
		if(o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0) return 1;

		else if(o1 == 0 && o3 != 0 && o4 != 0) return -1;
		else if(o2 == 0 && o4 != 0 && o3 != 0) return -2;
		else if(o3 == 0 && o1 != 0 && o2 != 0) return -3;
		else if(o4 == 0 && o2 != 0 && o1 != 0) return -4;
	}
	return 0; // doesn't include intersection
}

// calculates coordination at the intersection
Point calIntersection(Segment s1, Segment s2) {
	int flg = doIntersect(s1, s2);
	if(flg == 0) return {0, 0}; // no intersection
	else if(flg < 0) return {(double)flg, (double)flg}; // T-junction

	Point p1 = s1.left, q1 = s2.left, p2 = s1.right, q2 = s2.right;
	double S1 = ((q2.x-q1.x)*(p1.y-q1.y) - (q2.y-q1.y)*(p1.x-q1.x))/2.0;
	double S2 = ((q2.x-q1.x)*(q1.y-p2.y) - (q2.y-q1.y)*(q1.x-p2.x))/2.0;

	double x, y;
	x = p1.x + (p2.x-p1.x)*S1/(S1+S2);
	y = p1.y + (p2.y-p1.y)*S1/(S1+S2);

	return {x, y};
}

// generates Intersection
vector<Intersection> genInter(vector<Point> p, vector<Segment> arr) {
	set<Intersection> C;	
	vector<Intersection> retC;
	int M, c_M;

	M = arr.size();

	//FOR(j, i+1, M) {
	REP(i, M-1) {
		FOR(j, i+1, M) {
			Intersection tmpC;
			tmpC.p = calIntersection(arr[i], arr[j]);

			// detects intersection
			if(tmpC.p.x > 0) {
				auto itr = C.begin();
				for(; itr != C.end(); ++itr) {
					int cnt = 0;
					for(auto v : itr->s) if(i == v || j == v) ++cnt;
					if(cnt == 2) break;
				}
				if(itr == C.end()) {
					tmpC.s.insert(i);
					tmpC.s.insert(j);
					C.insert(tmpC);
				}
				// exists same intersection point
				else {
					itr->s.insert(i);
					itr->s.insert(j);
				}
			}
			else if(tmpC.p.x < 0) {
				auto itr = C.find(tmpC);
				if(itr != C.end()) {
					tmpC.s = itr->s;
					C.erase(itr);
				}
				if(tmpC.p.x == -1) {
					double pNum = (double)getIdx(p, arr[j].left);
					tmpC.p = {pNum, -1};
					tmpC.s.insert(i);
				}
				else if(tmpC.p.x == -2) {
					double pNum = (double)getIdx(p, arr[j].right);
					tmpC.p = {pNum, -1};
					tmpC.s.insert(i);
				}
				else if(tmpC.p.x == -3) {
					double pNum = (double)getIdx(p, arr[i].left);
					tmpC.p = {pNum, -1};
					tmpC.s.insert(j);
				}
				else if(tmpC.p.x == -4) {
					double pNum = (double)getIdx(p, arr[i].right);
					tmpC.p = {pNum, -1};
					tmpC.s.insert(j);
				}
				tmpC.flg = false;
				C.insert(tmpC);
			}
		}
	}
	for(auto itr : C) retC.push_back(itr);

	return retC;
}
// generates an adjacency matrix
// It can do execution with C included false
vector<vector<double>> genAdjMat(vector<Point> p, int N, vector<Segment> arr, vector<Intersection>& C) {	
	int M, c_num;
	vector<vector<double>> adj;
	unordered_map<int, vector<int>> mp;   // mp.first : segment number, mp.second : intersection numbers
					      // mp stores intersection numbers on the lines
	
	//N = p.size();
	M = arr.size();
	c_num = C.size();

	// It stores segment number on map(mp)
	int i=0;
	for(auto itr=C.begin(); itr!=C.end();) {
		if(itr->flg) {
			//p.push_back(itr->p);  // merge existed points and the intersection
		       	for(auto it=itr->s.begin(); it!=itr->s.end(); ++it) {
				//mp[*it].push_back(N+i);
				mp[*it].push_back(getIdx(p, itr->p));
			}
			++itr;
			++i;
		} else {
			for(auto it=itr->s.begin(); it!=itr->s.end(); ++it) {
				mp[*it].push_back((int)itr->p.x);
			}
			itr = C.erase(itr);
			--c_num;
		}
	}

	// if no intersection, push back NONE
	REP(i, arr.size()) {
		if(mp[i].empty()) mp[i].push_back(NONE);
	}

	// adj is initialize to -1
	REP(i, p.size()) {
		adj.emplace_back();
		REP(j, p.size()) adj[i].push_back(NONE);
	}

	// do ascending sort based on points to each vectors in map	
	REP(i, mp.size()) {
		REP(j, mp[i].size()-1) {
			FOR(k, j+1, mp[i].size()) {
				if(p[mp[i][j]] > p[mp[i][k]]) swap(mp[i][j], mp[i][k]);
			}
		}
	}

	for(auto itr=mp.begin(); itr!=mp.end(); ++itr) {
		auto it = itr->second.begin();

		// no intersection
		if(*it == NONE) {
			size_t idx_b = getIdx(p, arr[itr->first].left);
			size_t idx_e = getIdx(p, arr[itr->first].right);
			adj[idx_b][idx_e] = calDistance(arr[itr->first].left, arr[itr->first].right);
			adj[idx_e][idx_b] = calDistance(arr[itr->first].left, arr[itr->first].right);
			continue;
		}

		// calculates distance left point to 1st intersection
		size_t idx = getIdx(p, arr[itr->first].left);
		adj[idx][*it] = calDistance(arr[itr->first].left, p[*it]);
		adj[*it][idx] = calDistance(arr[itr->first].left, p[*it]);

		// calculates distance between intersections
		for(; it!=itr->second.end()-1; ++it) {
			adj[*it][*(it+1)] = calDistance(p[*it], p[*(it+1)]);
			adj[*(it+1)][*it] = calDistance(p[*it], p[*(it+1)]);
		}

		// calculates distance last intersection to right point
		idx = getIdx(p, arr[itr->first].right);
		adj[idx][*it] = calDistance(arr[itr->first].right, p[*it]);
		adj[*it][idx] = calDistance(arr[itr->first].right, p[*it]);
	}

	return adj;
}

Intersection findNearPoint(vector<Segment> arr, Point new_p) {
	double min = INF; // minimus distance
	int m;
	Intersection C;
	C.flg = false;

	REP(i, arr.size()) {
		// y = ax + b
		Point p, q, r;
		double a, b;
		double d;
		p = arr[i].left;
		q = arr[i].right;

		a = (p.y - q.y) / (p.x - q.x);
		b = p.y - (a * p.x);

		r.x = (a*(new_p.y - b) + new_p.x) / (a*a + 1);
		r.y = a*(a*(new_p.y - b) + new_p.x) / (a*a + 1) + b;

		if(r < p) r = p;
		else if(r > q) r = q;

		d = calDistance(new_p, r);

		if(d < min) {
			if(r == p || r == q) C.flg = false;
			else C.flg = true;
			min = d;
			C.p = r;
			m = i;
		}
	}
	C.s.insert(m);

	return C;
}

void addNewPoint(vector<Point>& p, vector<Segment>& arr, vector<Point> new_p, vector<Intersection>& C) {
	REP(i, new_p.size()) {
		Intersection tmpC;
		tmpC = findNearPoint(arr, new_p[i]);

		if(tmpC.p < new_p[i]) arr.push_back({tmpC.p, new_p[i]});
		else arr.push_back({new_p[i], tmpC.p});	

		auto c_itr = find(C.begin(), C.end(), tmpC);
		if(tmpC.flg) {
			// not exist the point
			if(c_itr == C.end()) {
				p.push_back(tmpC.p);
				C.push_back(tmpC);
			}
		}
	}
}

vector<int> getPath(vector<int> prev, int e) {
	vector<int> path;
	// from end to base point
	for(int now=e; now!=NONE; now=prev[now]) path.push_back(now);
	
	// reverse path bcz of little endian
	reverse(path.begin(), path.end());

	return path;
}

// adj : the adjacency matrix
//  b  : the base path
//  e  : the end point
PathInfo dijkstra(vector<vector<double>> adj, vector<int> b, int e) {
	int N = adj.size();
	vector<double> cost(N, INF);
	vector<bool> visited(N, false);
	// For restore the path
	vector<int> prev(N, NONE);
	// To get no visited node with min cost and
	priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> que;

	cost[b.front()] = 0;
	visited[b.front()] = true;
	for(auto itr = b.begin(); itr != b.end()-1; ++itr) {
		cost[*(itr+1)] = adj[*itr][*(itr+1)] + cost[*itr];
		visited[*(itr+1)] = true;
		prev[*(itr+1)] = *itr;
	}
	int now = b.back();

	while(1) {
		REP(i, adj[now].size()) {
			if(!visited[i] && adj[now][i] != NONE) {
				// if cost to the node decreases, cost is updated
				if(cost[i] > adj[now][i]+cost[now]) {
					cost[i] = adj[now][i] + cost[now];
					prev[i] = now;
				}
				que.push(make_pair(cost[i], i));
			}
		}
		if(que.empty()) break;
		now = que.top().second;
		que.pop();
		visited[now] = true;
	}

	vector<int> path;
	path = getPath(prev, e);


	return {cost[e], path};
}

vector<PathInfo> kth_ShortPath(vector<vector<double>> adj, int b, int e, int k) {
	vector<PathInfo> path_info; // determined paths (return values)
	set<PathInfo> tentPath; // tentative path
	
	path_info.push_back(dijkstra(adj, {b}, e)); // push back shortest path(1st)

	for(int i=2; i<=k; ++i) {
		PathInfo prevPath = path_info.back(); // (n-1)th path info when it calculates nth path

		vector<int> tmp_b; // searche each fragmented route in a row
		for(int j=0; j<prevPath.path.size()-1; ++j) {
			tmp_b.push_back(prevPath.path[j]);
			double tmp_dis;
			vector<pair<double, pair<int, int>>> er_edge; // erased edges
			
			// erase visited routes
			for(auto itr : path_info) {
				// check whether corresponding previous path to visited path
				if(checkHWords(itr.path, tmp_b)) {
					pair<int, int> route = make_pair(itr.path[j], itr.path[j+1]);

					if(adj[itr.path[j]][itr.path[j+1]] == NONE) continue;

					// store to the variable for restoring
					er_edge.push_back(make_pair(adj[itr.path[j]][itr.path[j+1]], route));
					// erase the road
					adj[itr.path[j]][itr.path[j+1]] = NONE;
					adj[itr.path[j+1]][itr.path[j]] = NONE;
				}
			}

			PathInfo dijPath = dijkstra(adj, tmp_b, e);
			if(dijPath.cost != INF) tentPath.insert(dijPath);

			// restore edges
			for(auto itr : er_edge) {
				adj[itr.second.first][itr.second.second] = itr.first;
				adj[itr.second.second][itr.second.first] = itr.first;
			}
		}
		// determine k-th shortest path
		if(tentPath.empty()) break;
		path_info.push_back(*tentPath.begin());
		tentPath.erase(tentPath.begin());
	}

	return path_info;
}

vector<vector<double>> adj;

// Trojan algorithm (lowlink) for finding bridges
vector<bool> visited;
vector<int> low, label;
int id = 0;

void dfs_init(int n) {
	id = 0;
	visited.clear();
	low.clear();
	label.clear();
	REP(i, n) {
		visited.push_back(false);
		low.push_back(NONE);
		label.push_back(NONE);
	}
}
void dfs(int at, int parent, vector<pair<int, int>>& bridge) {
	visited[at] = true;
	id++;
	low[at] = label[at] = id;
	
	REP(to, adj[at].size()) {
		// cannot visit or return to the original vertex
		if(adj[at][to] == NONE || to == parent) continue;
		if(!visited[to]) {
			dfs(to, at, bridge);
			low[at] = min(low[at], low[to]);
			if(label[at] < low[to]) bridge.push_back(make_pair(at, to));
		}
		else low[at] = min(low[at], low[to]);
	}
}

int main() {
	int N, M, P, Q;
	vector<Point> p;
	vector<Point> new_p;
	vector<Segment> arr;
	vector<string> b;
	vector<string> e;
	vector<int> k;
	//vector<vector<PathInfo>> path_info;

	cin >> N >> M >> P >> Q;
	REP(i, N) { 
		double x, y;
		cin >> x >> y;
		p.push_back({x, y});
	}

	REP(i, M) {
		int l, r;
		cin >> l >> r;
		arr.push_back({p[l-1], p[r-1]});
		if(arr[i].left>arr[i].right) swap(arr[i].left, arr[i].right);
	}

	REP(i, P) {
		double x, y;
		cin >> x >> y;
		new_p.push_back({x, y});
	}

	REP(i, Q) {
		string base, end;
		int k_th;
		cin >> base >> end >> k_th;

		b.push_back(base);
		e.push_back(end);
		k.push_back(k_th);
	}

	vector<Intersection> C;
	//set<Intersection> C;  // intersection
	
	int ordC_num = 0;
	C = genInter(p, arr);
	// the intersections are inserted into points
	for(auto itr : C) if(itr.flg) { p.push_back(itr.p); ++ordC_num; }

	for(auto itr : new_p) p.push_back(itr);
	addNewPoint(p, arr, new_p, C);

	//vector<vector<double>> adj;
	adj = genAdjMat(p, N, arr, C);

	cout << "-----intersection-------" << endl;
	if(ordC_num == 0) cout << "NA" << endl;
	else REP(i, ordC_num) {
		cout << p[N+i].x << " " << p[N+i].y << endl;
	}
	cout << "------------------------" << endl;

	cout << "-----new created point-----" << endl;
	REP(i, C.size()-ordC_num) {
		cout << p[N+ordC_num+P+i].x << " " << p[N+ordC_num+P+i].y << endl;
	}
	cout << "---------------------------" << endl;

	int c_num = C.size(); // number of intersection

	/*
	cout << "---------adjacency matrix----------" << endl;
	REP(i, adj.size()) {
		REP(j, adj[i].size()) cout << adj[i][j] << " ";
		cout << endl;
	}
	cout << "----------------------------------" << endl;*/

	// k-th shortest paths
	cout << "-----k-th shortest path-----" << endl;
	REP(i, Q) {
		vector<PathInfo> path_info;
		int bi, ei;
		bi = Point2Num(b[i], N, ordC_num);
		ei = Point2Num(e[i], N, ordC_num);
		if(bi > p.size()-ordC_num || ei > p.size()-ordC_num) {
			cout << "NA" << endl;
			continue;
		} else {
			path_info = kth_ShortPath(adj, bi, ei, k[i]);
			if(path_info[0].cost == INF) {
				cout << "NA" << endl;
				continue;
			}
			int s = path_info.size();
			REP(x, s) {
				cout << path_info[x].cost << endl;
				for(auto itr : path_info[x].path) cout << Num2Point(itr, N, ordC_num, P) << " ";
				cout << endl;
			}
		}
	}
	cout << "----------------------------" << endl;

	// print all points
	//for(auto itr : p) cout << Num2Point(getIdx(p, itr), N, ordC_num, P) << " : " << itr.x << " " << itr.y << endl;
	
	// calculates bridges in the graph
	vector<pair<int, int>> bridge;
	dfs_init(adj.size());
	cout << "-----bridge-----" << endl;
	REP(i, p.size()) {
		vector<pair<int, int>> bridge;	
		if(!visited[i]) dfs(i, i, bridge);
		for(auto itr : bridge) cout << Num2Point(itr.first, N, ordC_num, P) << " " << Num2Point(itr.second, N, ordC_num, P) << endl;
	}	
	//cout << "-----bridge-----" << endl;
	//for(auto itr : bridge) cout << Num2Point(itr.first, N, ordC_num, P) << " " << Num2Point(itr.second, N, ordC_num, P) << endl;
	cout << "----------------" << endl;
}
