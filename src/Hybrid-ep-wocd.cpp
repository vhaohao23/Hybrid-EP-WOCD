#include<bits/stdc++.h>
using namespace std;
int pop=100;
const int T=100;
int N;
int NE;
const double p=1.0;
const double lenP=5.0;
int macom=0;
vector<vector<bool>> A;
vector<vector<int>> x(pop+1);
vector<vector<int>> e;
vector<int> xBest;
vector<int> d;
vector<vector<int>> dk(pop+1);
vector<vector<int>> lk(pop+1);

random_device rd;   
mt19937 gen(rd()); 

vector<queue<int>> remainComs(pop+1);

const double eliminateRate=0.4;

double modularity(vector<int> dk,vector<int> lk){
    double Q=0;

    for (int i=1;i<=N;i++){
        Q+=double(lk[i])/double(NE)-pow(double(dk[i])/double(2*NE),2.0);
    }

    return Q;
}

void transfer(vector<int> &dk,vector<int> &lk,vector<int> l,int i,int l1,int l2){
    dk[l1]-=d[i];
    dk[l2]+=d[i];
    for (int v:e[i]){
        if (l1==l[v])
            --lk[l1];
        if (l2==l[v])
            ++lk[l2];
    }
}

void LAR_rand(vector<vector<int>> &a){
    for (int u=1;u<=N;u++){
        if (!e[u].size()) continue;
        uniform_int_distribution<int> disv(0,e[u].size()-1);
        int v=e[u][disv(gen)];
        a[u].push_back(v);
        a[v].push_back(u);
    }
}

bool Excute(int &ll1,int &lf1,int &ll2,int &lf2,int &dl1,int &df1,int &dl2,int &df2){
    double deltaQ=double(ll1-lf1+ll2-lf2)/double(NE)-(pow(double(dl1-df1),2.0)+pow(double(dl2-df2),2.0))/double(4*NE*NE);
    
    if (deltaQ<=0){
        ll1=lf1;
        ll2=lf2;
        dl1=df1;
        dl2=df2;
        return 0;
    }
    return 1;
}

vector<int> decoding(vector<vector<int>> a){
    bool dd[N+1]={};
    vector<int> l(N+1);
    int cnt=0;

    for (int i=1;i<=N;i++)
        if (!dd[i]){
            ++cnt;
            queue<int> q;
            q.push(i);
            while (!q.empty()){
                int u=q.front();
                q.pop();
                l[u]=cnt;
                for (int v:a[u])
                    if (!dd[v]){
                        dd[v]=true;
                        q.push(v);
                    }
            }
        }
    return l;
}

void initialization(){
    int s;
    for (int p=1;p<=pop;p++){
        vector<vector<int>> a(N+1);
        LAR_rand(a);
        x[p]=decoding(a);         
        
        s=*max_element(x[p].begin(),x[p].end());
        dk[p].resize(N+1,0);
        lk[p].resize(N+1,0);                            
        
        for (int u=1;u<=N;u++){
            dk[p][x[p][u]]+=d[u];
             
            for (int v:e[u])
                if (x[p][u]==x[p][v]&&u<v){
                    ++lk[p][x[p][u]];
                }
        }

        for (int i=s+1;i<=N;i++)
            remainComs[p].push(i);
    }    
}

void movingToPrey(vector<int> &l,vector<int> &dk,vector<int> &lk,double k,queue<int> &remainCom){
    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++){
        transfer(dk,lk,l,pos[i],l[pos[i]],xBest[pos[i]]);
        l[pos[i]]=xBest[pos[i]];
    }
}

void randomWalk(vector<int> &l,vector<int> &dk,vector<int> &lk,double k,queue<int> &remainCom){
    vector<int> ranPop;
    for (int i=1;i<=pop;i++)
        ranPop.push_back(i);
    shuffle(ranPop.begin(),ranPop.end(),gen);   

    vector<vector<int>> P;
    for (int i=0;i<=lenP-1;i++)
        P.push_back(x[ranPop[i]]);
    

    vector<int> randNode;
    for (int i=1;i<=N;i++)
        randNode.push_back(i);
    shuffle(randNode.begin(),randNode.end(),gen);


    for (int i=lenP-1;i>=0;i--){
        int kCommma=int(k)/(i+1)+int(int(k)%(i+1)>0);

        while (kCommma--){
            --k;
            transfer(dk,lk,l,randNode[k],l[randNode[k]],P[i][randNode[k]]);
            l[randNode[k]]=P[i][randNode[k]];
        }
    }
}

void encirlingThePrey(vector<int>&l,vector<int> &dk,vector<int> &lk,double r,queue<int> &remainCom){
    int k=int(r*double(N));

    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++){
        transfer(dk,lk,l,pos[i],l[pos[i]],xBest[pos[i]]);
        l[pos[i]]=xBest[pos[i]];
    }
}

void mutation(vector<int> &l,vector<int> &dk,vector<int> &lk,double u,queue<int> &remainCom){
    uniform_real_distribution<double> dis(0,1);
    
    int S = *max_element(l.begin(), l.end());
    vector<int> ltmp;
    vector<int> dktmp;
    vector<int> lktmp;
    for (int i=1;i<=N;i++){
        double x=dis(gen);
        if (x<u){
            ltmp=l; 
            dktmp=dk;
            lktmp=lk;

            S=remainCom.front();
            remainCom.pop();

            double y=dis(gen);
            if (y<0.5){
                transfer(dk,lk,l,i,l[i],S);                
                l[i]=S;
                
            }
            else{
                transfer(dk,lk,l,i,l[i],S);
                l[i]=S;
                

                for (int neigbor:e[i]){
                    transfer(dk,lk,l,neigbor,l[neigbor],S);
                    l[neigbor]=S;
                }
            }
            
            if (modularity(dk,lk)<modularity(dktmp,lktmp)){
                remainCom.push(S);
                l=ltmp;
                dk=dktmp;
                lk=lktmp;
            }
            macom=max(macom,S);
        }
    }
}

bool isBoundaryNode(vector<int> l,int u){
    for (int v:e[u])
        if (l[u]!=l[v])
            return true;
    return false;
}

void boudaryNodeAdjustment(vector<int> &l,vector<int> &dk,vector<int> &lk,queue<int> &remainCom){
    vector<int> tmpl;
    vector<int> dktmp;
    vector<int> lktmp;
    int s=*max_element(l.begin(),l.end());
    bool dd[s+1]={};
    
    for (int i=1;i<=N;i++){
        if (isBoundaryNode(l,i)){
            dd[l[i]]=true;
            for (int neighbor:e[i])
                if (l[i]!=l[neighbor] && !dd[l[neighbor]]){
                    dd[l[neighbor]]=true;
                    
                    tmpl=l;
                    dktmp=dk;
                    lktmp=lk;

                    transfer(dk,lk,l,i,l[i],l[neighbor]);
                    l[i]=l[neighbor];
                    
                    if (modularity(dk,lk)<modularity(dktmp,lktmp)){
                        dk=dktmp;
                        lk=lktmp;
                        l=tmpl;
                    }

                    // if (!Excute(lk[tmpl[i]],lktmp[tmpl[i]],lk[l[neighbor]],lktmp[l[neighbor]],dk[tmpl[i]],dktmp[tmpl[i]],lk[l[neighbor]],lktmp[l[neighbor]] ) )
                    //     l=tmpl;
                }

            for (int i=1;i<=s;i++)
                dd[i]=false;
        }
    }
}

void mimic(vector<int> &lLose,vector<int>&dkLose,vector<int> &lkLose, vector<int> &lWin,vector<int>&dkWin,vector<int> &lkWin){
    uniform_int_distribution dis(1,N);
    int v=dis(gen);

    for (int i=1;i<=N;i++)
        if (lLose[i]==lLose[v]){
            transfer(dkLose,lkLose,lLose,i,lLose[i],xBest[i]);
            lLose[i]=xBest[i];
        }
}

void EPD(){
    if (x.size()<10) return;

    vector<pair<double, int>> modularityValues;
    for (int i = 1; i <= pop; i++) {
        double modValue = modularity(dk[i],lk[i]);  
        modularityValues.push_back({modValue, i});
    }

    sort(modularityValues.begin(), modularityValues.end());
    reverse(modularityValues.begin(),modularityValues.end());

    vector<vector<int>> sortedX(pop + 1);
    vector<vector<int>> sorteddk(pop + 1);
    vector<vector<int>> sortedlk(pop + 1); 
    for (int i = 0; i < pop; i++) {
        sortedX[i + 1] = x[modularityValues[i].second];
        sorteddk[i + 1] = dk[modularityValues[i].second];
        sortedlk[i + 1] = lk[modularityValues[i].second];
    }

    x=sortedX;
    dk=sorteddk;
    lk=sortedlk;

    int numCompetitor=pow(2, int( log2(pop) ) );
   
    // for (int p=pop;p>=numCompetitor+1;p--){
    //     x.erase(x.begin() + p);
    //     dk.erase(dk.begin() + p);
    //     lk.erase(lk.begin() + p);   
    //     --pop;
    // }

    //tourament
    vector<int> idx(pop);
    iota(idx.begin(),idx.end(),1);
    shuffle(idx.begin(), idx.end(), gen);
    vector<vector<int>> newx(pop + 1);
    vector<vector<int>> newdk(pop + 1);
    vector<vector<int>> newlk(pop + 1);
    for (int i = 0; i < pop; i++) {
        newx[i + 1] = x[idx[i]];
        newdk[i + 1] = dk[idx[i]];
        newlk[i + 1] = lk[idx[i]];
    }
    x = newx;
    dk = newdk;
    lk = newlk;
    uniform_real_distribution disEli(0.0,1.0);

    for (int p=1;p<pop;p+=2){
        // cout<<p<<" "<<pop<<"\n";
        if (modularity(dk[p],lk[p]) > modularity(dk[p+1],lk[p+1])){
            // cout<<1<<"\n";
            if (disEli(gen)>eliminateRate)
                mimic(x[p+1],dk[p+1],lk[p+1],x[p],dk[p],lk[p]);
            else{
                x.erase(x.begin() + p+1);
                dk.erase(dk.begin() + p+1);
                lk.erase(lk.begin() + p+1);
                --pop;
                --p;
            }
        }
        else {
            // cout<<1<<"\n";
            if (disEli(gen)>eliminateRate)
                mimic(x[p],dk[p],lk[p],x[p+1],dk[p+1],lk[p+1]);
            else {
                x.erase(x.begin() + p);
                dk.erase(dk.begin() + p);
                lk.erase(lk.begin() + p);
                --pop;
                --p;
            }
        }
    }
}

void updateLocation(vector<int> &l,int t,vector<int> &dk,vector<int> &lk,queue<int> &remainCom){
    uniform_real_distribution<double> dis(0,1);
    double alpha=dis(gen),beta=dis(gen);
    if (alpha<0.5){
        double k=p*double(t)*double(N)/double(T);

        if (beta<0.5)
            movingToPrey(l,dk,lk,k,remainCom);
        else
            randomWalk(l,dk,lk,k,remainCom);
    }
    else {
        uniform_real_distribution<double> disl(-1,1);
        double ll=disl(gen);
        double r=abs(cos(2*M_PI*ll));
        encirlingThePrey(l,dk,lk,r,remainCom);
    }
}

void EP_WOCD(){
    initialization();
    double ans=0;
    for (int i=1;i<=pop;i++)
        if (modularity(dk[i],lk[i])>ans){
            ans=modularity(dk[i],lk[i]);
            xBest=x[i];
        }

    bool dd[N+1]={};
    for (int t=1;t<=T;t++){
        for (int p=1;p<=pop;p++){
            updateLocation(x[p],t,dk[p],lk[p],remainComs[p]);
            mutation(x[p],dk[p],lk[p],0.3,remainComs[p]);
            boudaryNodeAdjustment(x[p],dk[p],lk[p],remainComs[p]);
            
        }
        for (int i=1;i<=pop;i++){
            if (modularity(dk[i],lk[i])>ans){
                ans=modularity(dk[i],lk[i]);
                xBest=x[i];
            }
        }

        for (int p=1;p<=pop;p++){
            while (remainComs[p].size())
            {
                remainComs[p].pop();
            }

            for (int i=1;i<=N;i++) dd[i]=0;
            for (int i=1;i<=N;i++) dd[x[p][i]]=1;
            for (int i=1;i<=N;i++)
                if (!dd[i]) remainComs[p].push(i); 
        }

        
        EPD();   
        // cout<<pop<<"\n"  ;
    }    

    cout<<ans<<"\n";
    for (int i=1;i<=N;i++)
        cout<<xBest[i]<<" ";
}
int main(){
    clock_t tStart = clock();

    freopen("input.txt","r",stdin);

    cin>>N;
    cin>>NE;

    d.resize(N+1);
    A.resize(N + 1, vector<bool>(N + 1, 0));
    e.resize(N+1);
    for (int i=1;i<=NE;i++){
        int u,v;
        cin>>u>>v;
        e[u].push_back(v);
        e[v].push_back(u);
        d[u]++,d[v]++;
        A[u][v]=A[v][u]=true;
    }

    EP_WOCD();

    printf("\nTime taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}
