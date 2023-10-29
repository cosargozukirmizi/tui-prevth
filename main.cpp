#include <iostream>
#include <stdlib.h>  // for EXIT_SUCCESS
#include <functional>
#include <memory>    // for allocator, __shared_ptr_access
#include <string>  // for string, operator+, basic_string, to_string, char_traits
#include <vector>  // for vector, __alloc_traits<>::value_type
#include "ftxui/component/captured_mouse.hpp"  // for ftxui
#include "ftxui/component/component.hpp"  // for Menu, Renderer, Horizontal, Vertical
#include "ftxui/component/component_base.hpp"  // for ComponentBase
#include "ftxui/component/screen_interactive.hpp"  // for Component, ScreenInteractive
#include "ftxui/dom/elements.hpp"  // for text, Element, operator|, window, flex, vbox
#include <map>
#include <gmpxx.h>
#include <iomanip>
#include <thread>
#include <semaphore>
#include <algorithm>
#include <cassert>

using namespace ftxui;


Component Window(std::string title, Component component) {
  return Renderer(component, [component, title] {  //
    return window(text(title), component->Render());
  });
}


Component Window2(std::string title, Component component) {
  return Renderer(component, [component, title] {  //
    return window(text(title), component->Render()) |vscroll_indicator|frame |
           size(HEIGHT, LESS_THAN, 10) | border;
  });
}

ButtonOption ButtonStyle() {
  auto option = ButtonOption::Animated();
  option.transform = [](const EntryState& s) {
    auto element = text(s.label);
    if (s.focused) {
      element |= bold;
    }
    return element | center | borderEmpty | flex;
  };
  return option;
}


void Settings(std::vector<std::vector<mpq_class>>& vd, std::vector<std::vector<std::string>> menu_entries, int menu_selected[]) {
  auto screen = ScreenInteractive::TerminalOutput();
  auto back_button = Button("Back", screen.ExitLoopClosure());
  std::string odename = menu_entries[0][menu_selected[0]];
  auto goto_1 = Button("Add 1/2 to " + odename + " f0 initVal", [&vd, menu_selected] { (vd.at(menu_selected[0])).at(0)=(vd.at(menu_selected[0])).at(0)+mpq_class(1,2); });
  auto goto_2 = Button("Subtract 1/2 from " + odename + " f0 initVal", [&vd, menu_selected] { (vd.at(menu_selected[0])).at(0)=(vd.at(menu_selected[0])).at(0)-mpq_class(1,2); });
  auto goto_3 = Button("Add 1/2 to " + odename + " f1 initVal", [&vd, menu_selected] { (vd.at(menu_selected[0])).at(1)=(vd.at(menu_selected[0])).at(1)+mpq_class(1,2); });
  auto goto_4 = Button("Subtract 1/2 from "+ odename + " f1 initVal", [&vd, menu_selected] { (vd.at(menu_selected[0])).at(1)=(vd.at(menu_selected[0])).at(1)-mpq_class(1,2); });
  auto layout = Container::Vertical({
      back_button,
      goto_1,
      goto_2,
      goto_3,
      goto_4,
  });
  auto renderer = Renderer(layout, [&] {
    return vbox({
//               text("path: " + path),
               separator(),
               back_button->Render(),
               goto_1->Render(),
               goto_2->Render(),
               goto_3->Render(),
               goto_4->Render(),
           }) |
           border;
  });
  screen.Loop(renderer);
}


void extendSpace (std::string& xx, const std::vector<int>& equationVector, std::vector<mpq_class>& myCoeffs, std::vector<mpq_class>& initVal, const int max_iter, const int f_num, const int p_prec, const int t_maxChoice);
void construct2F(std::vector<int>& rowInd, std::vector<int>& colInd, std::vector<mpq_class>& myValue, const std::vector<int>& runVec, const std::map<std::vector<int>, int>& myMap, const std::vector<mpq_class>& myCoeffs, const int stEq);
void condensedKroneckerProduct(std::vector<mpq_class>& result, const std::vector<mpq_class>& a, const std::vector<mpq_class>& b);
void sparseMatTimesVec(std::vector<mpq_class>& result, const std::vector<int>& rowInd, const std::vector<int>& colInd, const std::vector<mpq_class>& myCoeffs, const std::vector<mpq_class>& x);
void constructAugmentedInitVal(std::vector<mpq_class>& runInitVal, std::vector<mpq_class>& initVal, const std::map<std::vector<int>, int>& myMap);


int main() {


  auto screen = ScreenInteractive::TerminalOutput();
  screen.TrackMouse(false);

  const int max_iter = 3;

  std::vector<std::vector<int>> allODEs;

  std::vector<std::vector<mpq_class>> allCoeffs;

  std::vector<std::vector<mpq_class>> allInitVals;


  const std::vector<int> vanderPol
  {
    1, 0,  0, 0,  1, 0,
    1, 0,  1, 0,  2, 0,
    1, 0,  0, 0,  0, 1,
    0, 1,  0, 0,  1, 0
  };

  allODEs.push_back(vanderPol);

  std::vector<mpq_class> vanderPolCoeffs{
   mpq_class(1,1), mpq_class(-1,3), mpq_class(-1,1), mpq_class(1,1)
  };

  allCoeffs.push_back(vanderPolCoeffs);

  std::vector<mpq_class> vanderPolInitVal{
   mpq_class(1,2), mpq_class(1,2)
  };

  allInitVals.push_back(vanderPolInitVal);


  const std::vector<int> quarticAnharmonicOscillator
  {
    1, 0,  0, 0,  0, 1,
    0, 1,  0, 0,  1, 0,
    0, 1,  1, 0,  2, 0
  };

  allODEs.push_back(quarticAnharmonicOscillator);


  std::vector<mpq_class> quarticAnharmonicOscillatorCoeffs{
   mpq_class(1), mpq_class(-1), mpq_class(-1)
  };   // mu is 1, k1 is 1, k2 is 1

  allCoeffs.push_back(quarticAnharmonicOscillatorCoeffs);


  const std::vector<mpq_class> quarticAnharmonicOscillatorInitVal{
   mpq_class(1,2),mpq_class(1,2)
  };

  allInitVals.push_back(quarticAnharmonicOscillatorInitVal);

  const std::vector<int> henonHeiles
  {
    1, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,
    0, 1, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0,
    0, 1, 0, 0,  0, 0, 1, 0,  1, 0, 0, 0,
    0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 1,
    0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 1, 0,
    0, 0, 0, 1,  1, 0, 0, 0,  1, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 1, 0,  0, 0, 1, 0
  };

  allODEs.push_back(henonHeiles);


  std::vector<mpq_class> henonHeilesCoeffs{
   mpq_class(1), mpq_class(-1), mpq_class(-2), mpq_class(1), mpq_class(-1), mpq_class(-1), mpq_class(1)
  }; // lambda is 1

  allCoeffs.push_back(henonHeilesCoeffs);


  const std::vector<mpq_class> henonHeilesInitVal{
   mpq_class(1,2), mpq_class(1,2), mpq_class(1,2), mpq_class(1,2)
  };

  allInitVals.push_back(henonHeilesInitVal);

  const std::vector<int> rabinovichFabrikant
  {
    1, 0, 0,  0, 0, 1,  0, 1, 0,
    1, 0, 0,  0, 0, 0,  0, 1, 0,
    1, 0, 0,  1, 0, 0,  1, 1, 0,
    1, 0, 0,  0, 0, 0,  1, 0, 0,
    0, 1, 0,  0, 0, 1,  1, 0, 0,
    0, 1, 0,  0, 0, 0,  1, 0, 0,
    0, 1, 0,  1, 0, 0,  2, 0, 0,
    0, 1, 0,  0, 0, 0,  0, 1, 0,
    0, 0, 1,  0, 0, 0,  0, 0, 1,
    0, 0, 1,  0, 1, 0,  1, 0, 1
  };

  allODEs.push_back(rabinovichFabrikant);


  std::vector<mpq_class> rabinovichFabrikantCoeffs{
   mpq_class(1), mpq_class(-1), mpq_class(1), mpq_class(1), mpq_class(3), mpq_class(1), mpq_class(-1), mpq_class(1), mpq_class(-2), mpq_class(-2)
  }; // gamma is 1, alpha is 1


  allCoeffs.push_back(rabinovichFabrikantCoeffs);

  const std::vector<mpq_class> rabinovichFabrikantInitVal{
   mpq_class(1,2), mpq_class(1,2), mpq_class(1,2)
  };

  allInitVals.push_back(rabinovichFabrikantInitVal);



  int menu_selected[] = {0, 0, 0, 0};
  std::vector<std::vector<std::string>> menu_entries = {
      {
          "vdPol",
          "QuartAnhOsc",
          "HenHeil",
          "RabFab",
      },
      {
          "1/10",
          "1/8",
          "1/6",
          "1/4",
          "1/2",
          "1",
      },
      {
          "f0",
          "f1",
      },
      {
          "100",
          "120",
          "140",
      }
  };

  int num_iter = 3;


  auto button_settings = Button("Settings", [&allInitVals, &menu_entries, &menu_selected] { Settings(allInitVals, menu_entries, menu_selected); });
  auto button_p1 = Button("num_iter+1", [&num_iter] { num_iter = num_iter+1; });
  auto button_m1 = Button("num_iter-1", [&num_iter] { if(num_iter>3) {num_iter = num_iter - 1;} });
  auto button_p10 = Button("num_iter+10", [&num_iter] { num_iter = num_iter+10; });
  auto button_m10 = Button("num_iter-10", [&num_iter] { if(num_iter>12) {num_iter = num_iter-10;} });
  auto button_quit = Button("Quit", screen.ExitLoopClosure());

  auto menu_global = Container::Vertical(
      {
          Window("ODE preset", Radiobox(&menu_entries[0], &menu_selected[0])),
          Window("t_max", Radiobox(&menu_entries[1], &menu_selected[1])),
          button_p1,
          button_m1,
          button_p10,
          button_m10,
          Window("f_num", Radiobox(&menu_entries[2], &menu_selected[2])),
          Window("p_prec", Radiobox(&menu_entries[3], &menu_selected[3])),
          button_settings,
          button_quit,
      });

  auto info = Renderer([&] {
    std::string xx = "";
    extendSpace (xx, allODEs[menu_selected[0]], allCoeffs[menu_selected[0]], allInitVals[menu_selected[0]], num_iter, menu_selected[2], menu_selected[3], menu_selected[1]);

    return window(text("Output"),
                  vbox({
                      text("num_iter	= " +
                           std::to_string(num_iter)),
                      paragraph("\n"),
                      paragraph(xx),
                  })) | flex;
  });


  auto global = Container::Horizontal({
      menu_global,
      info,
  });


  screen.Loop(global);
  return EXIT_SUCCESS;
}


void extendSpace (std::string& xx, const std::vector<int>& equationVector, std::vector<mpq_class>& myCoeffs, std::vector<mpq_class>& initVal, const int max_iter, const int f_num, const int p_prec, const int t_maxChoice)    // time value under consideration
{
  using namespace std;

  int p_precAsNum = 100;

  mpq_class t_max(1,10);

  if ( t_maxChoice == 1)
    t_max = (t_max * 10) / 8;
  else if ( t_maxChoice == 2)
    t_max = (t_max * 10) /6; 
  else if ( t_maxChoice == 3)
    t_max = (t_max * 10) / 4; 
  else if ( t_maxChoice == 4)
    t_max = (t_max * 10) / 2; 
  else if ( t_maxChoice == 5)
    t_max = t_max  * 10; 


  if( p_prec == 1 )
  {
    p_precAsNum = 120;
  }
  else if( p_prec == 2 )
  {
    p_precAsNum = 140;
  }


  assert(max_iter > 1);
  const int stEq = initVal.size();
  vector<int> runVec = equationVector;
  map<vector<int>, int> myMap;

  bool appears = 0;

  vector<int>::size_type i = 0;
  auto numEqs = stEq;

  vector<int> rightHandSide (stEq, 0);
  vector<int> tempLeft (stEq, 0);

  for (vector<int>::size_type i = 0; i < runVec.size (); i+=3*stEq)
  {
    tempLeft.assign(equationVector.begin()+i, equationVector.begin()+i+stEq);
    myMap.insert({tempLeft, 0});
  }

  while (i < runVec.size ())
  {
    for (auto j = stEq; j <= 2 * stEq; j += stEq)
    {
      for (auto k = 0; k < stEq; k++)
      {
        rightHandSide[k] = runVec[i + j + k];
      }
      if (myMap.find(rightHandSide) == myMap.end())
      {
         myMap.insert({rightHandSide, 0});
      }

      for (vector<int>::size_type i = 0; i < runVec.size (); i += 3 * stEq)
      {
        vector<int> subvector(stEq, 0);
        subvector = {runVec.begin () + i, runVec.begin () + i + stEq};
        if (subvector == rightHandSide)
        {
          appears = 1;
          break;
        }
      }

      if (appears == 0)
      {
        ++numEqs;
        vector<int> allZero(stEq,0);
        vector<int> zeroEntry(3*stEq, 0);

        if (rightHandSide == allZero)
        {
          runVec.insert (runVec.end (), zeroEntry.begin(), zeroEntry.end() );
          myCoeffs.push_back(mpq_class(0));
        }
        else
        {
          for (vector<int>::size_type i = 0; i < equationVector.size (); i += 3*stEq)
          {
            vector<int>::const_iterator it;
            it = find (equationVector.begin () + i, equationVector.begin () + i + 3*stEq, 1);
            vector<int>::size_type oneInd = (it - equationVector.begin ());

            bool x = 0;
            vector<int> myTemp(stEq, 0);

            for (auto j = 0; j < stEq; j++)
            {
              if (oneInd == i+j && rightHandSide[j] > 0)
              {
                x = 1;
                for (auto k = 0 ; k < stEq; k++)
                {
                  myTemp[k] =
                    equationVector[i + stEq + k] + equationVector[i + stEq + stEq + k] +
                    rightHandSide[k];
                  if (j == k)
                  {
                    myCoeffs.push_back(rightHandSide[j] * myCoeffs[i/(3*stEq)]);
                    myTemp[k] -= 1;
                  }
                }
              }
            }

            if (x == 0)
            {
              continue;
            }
            vector<int> lefts;
            int oneCount = 0;

            for (auto counter = 0; counter < stEq; counter++)
              runVec.push_back (rightHandSide[counter]);

            for (vector<int>::size_type i = 0; i < myTemp.size (); ++i)
            {
              if (myTemp[i] == 1)
              {
                ++oneCount;
              }
              if (myTemp[i] != 1 || (myTemp[i]==1 && (oneCount%2==1)))
              {
                lefts.push_back(myTemp[i] / 2);
                runVec.push_back(myTemp[i] / 2);
              }
              else if (myTemp[i]==1 && oneCount%2==0)
              {
                lefts.push_back(1);
                runVec.push_back(1);
              }
            }
            for (vector<int>::size_type i = 0; i < lefts.size (); ++i)
            {
              runVec.push_back(myTemp[i] - lefts[i]);
            }
          }
        }
      }
      appears = 0;
    }

    i += 3*stEq;
  }


  for(auto myDummyVariable = 0; auto& [key, value] : myMap)
  {
    value = myDummyVariable++;
  }


  vector<mpq_class> runInitVal(myMap.size());


  constructAugmentedInitVal(runInitVal, initVal, myMap);


  vector<int> rowInd;
  vector<int> colInd;
  vector<mpq_class> myValue;

  construct2F(rowInd, colInd, myValue, runVec, myMap, myCoeffs, stEq);

  vector<vector<mpq_class>> rho(max_iter, vector<mpq_class>(runInitVal.size()));

  rho.at(0) = (runInitVal);

  vector<mpq_class> resOfKronProd(runInitVal.size()*(runInitVal.size()+1)/2, 0);

  condensedKroneckerProduct(resOfKronProd, rho[0], rho[0]);
  std::transform(resOfKronProd.begin(), resOfKronProd.end(), resOfKronProd.begin(), [](mpq_class& c){return c / 2;});

  sparseMatTimesVec(rho.at(1), rowInd, colInd, myValue, resOfKronProd);

  std::binary_semaphore smphSignalMainToThread{1}, smphSignalThreadToMain{0};  // oddSummer has green light, it is good to go

  auto evenSummer = [&]()
  {

    for(auto j=2; j < max_iter-1; j=j+2)
    {
       vector<mpq_class> resOfKronProd(runInitVal.size()*(runInitVal.size()+1)/2, 0);
       vector<mpq_class> resOfMatVecProd(runInitVal.size(), 0);
       vector<mpq_class> resOfKronProdAccumulator(runInitVal.size()*(runInitVal.size()+1)/2, 0);

       std::fill(resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), 0);
       condensedKroneckerProduct(resOfKronProd, rho[j/2], rho[j/2]);
       std::transform(resOfKronProd.begin(), resOfKronProd.end(), resOfKronProd.begin(), [](mpq_class& c){return c / 2;});
       std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());


       for(auto k=(j/2-1); k>0; --k)
       {
         condensedKroneckerProduct(resOfKronProd, rho[k], rho[j-k]);
         std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());
       }

       smphSignalThreadToMain.acquire();

       condensedKroneckerProduct(resOfKronProd, rho[0], rho[j]);
       std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());

       sparseMatTimesVec(resOfMatVecProd, rowInd, colInd, myValue, resOfKronProdAccumulator);

       mpq_class myConst(j+1);

       std::transform(resOfMatVecProd.begin(), resOfMatVecProd.end(), resOfMatVecProd.begin(), [&myConst](mpq_class& c){return c / myConst;});

       rho.at(j+1) = resOfMatVecProd;

       smphSignalMainToThread.release();
    }

  };


  auto oddSummer = [&]()
  {
    for(auto j=1; j < max_iter-1; j=j+2)
    {
       vector<mpq_class> resOfKronProd(runInitVal.size()*(runInitVal.size()+1)/2, 0);
       vector<mpq_class> resOfMatVecProd(runInitVal.size(), 0);
       vector<mpq_class> resOfKronProdAccumulator(runInitVal.size()*(runInitVal.size()+1)/2, 0);

       std::fill(resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), 0);

       for(auto k=(j-1)/2; k>0; --k)  // 0 not included for a reason
       {
         condensedKroneckerProduct(resOfKronProd, rho[k], rho[j-k]);
         std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());
       }

       smphSignalMainToThread.acquire();

       condensedKroneckerProduct(resOfKronProd, rho[0], rho[j]);
       std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());

       sparseMatTimesVec(resOfMatVecProd, rowInd, colInd, myValue, resOfKronProdAccumulator);

       mpq_class myConst(j+1);

       std::transform(resOfMatVecProd.begin(), resOfMatVecProd.end(), resOfMatVecProd.begin(), [&myConst](mpq_class& c){return c / myConst;});

       rho.at(j+1) = resOfMatVecProd;

       smphSignalThreadToMain.release();
    }
  };

  thread t(oddSummer);
  evenSummer();

  t.join();


  const mpq_class t_step(t_max/10);


  for(mpq_class indQ=0; indQ<=t_max; indQ=indQ+t_step)
  {

    vector<mpq_class> solution(runInitVal.size(), 0);

  //  cout << '\n';
    mpq_class leftPart(1);

    for(auto i = 0; i < max_iter; ++i)
    {
      vector<mpq_class> solutionIngredient(rho[i].begin(), rho[i].end());

      std::transform (solutionIngredient.begin(), solutionIngredient.end(), solutionIngredient.begin(), [&leftPart](auto& c){return c*leftPart;});
      std::transform (solution.begin(), solution.end(), solutionIngredient.begin(), solution.begin(), std::plus<mpq_class>());


      // truncation order is the number of terms.
      // truncation order is also the number of iterations.
      if( i == max_iter - 1)
      {

        vector<int> myKey (stEq, 0);
        myKey.at(f_num) = 1;

        auto myValue = myMap.at(myKey);

        char tempStr[p_precAsNum*10];
        string tempStr2;

        mpf_class f(solution[myValue], p_precAsNum*10);
        gmp_sprintf (tempStr, "%.*Ff\n", p_precAsNum, f);
        tempStr2 = tempStr;
        xx = xx + tempStr2;
        xx = xx + " ";
      }
      leftPart *= indQ;
    }

  }
  xx = xx + " ";

  return;
}


int findPlace(const int a, const int b, const int n)
{
  return a*n+b-(a*(a+1))/2;
}


void constructAugmentedInitVal(std::vector<mpq_class>& runInitVal, std::vector<mpq_class>& initVal, const std::map<std::vector<int>, int>& myMap)
{
  using namespace std;

  for(const auto& [key, value] : myMap)
  {
    mpq_class tempVal(1);

    for(vector<int>::size_type i = 0; i<key.size(); ++i)
    {
      mpq_class tempInside(1);

      for(auto j = 1; j<=key[i]; ++j)
      {
        tempInside = tempInside * initVal[i];
      }
      tempVal = tempVal * tempInside;
    }

    runInitVal[value] = tempVal;
  }
}


// The function below constructs 2*F.
void construct2F(std::vector<int>& rowInd, std::vector<int>& colInd, std::vector<mpq_class>& myValue, const std::vector<int>& runVec, const std::map<std::vector<int>, int>& myMap, const std::vector<mpq_class>& myCoeffs, const int stEq)
{
  using namespace std;

  vector<int> left(stEq);
  vector<int> middle(stEq);
  vector<int> right(stEq);

  for(vector<int>::size_type i = 0; i < runVec.size(); i+=3*stEq)
  {
    left = { runVec.begin() + i, runVec.begin() + i + stEq };
    middle = { runVec.begin() + i + stEq, runVec.begin() + i + stEq +stEq };
    right = { runVec.begin() + i + stEq + stEq, runVec.begin() + i + stEq + stEq + stEq};

    const int leftAsInd = myMap.at(left);
    const int middleAsInd = myMap.at(middle);
    const int rightAsInd = myMap.at(right);

    rowInd.push_back(leftAsInd);

    colInd.push_back(findPlace(middleAsInd, rightAsInd, myMap.size()));
    if (middleAsInd == rightAsInd)
    {
      myValue.push_back(myCoeffs[i/(3*stEq)]*2);
    }
    else
    {
      myValue.push_back(myCoeffs[i/(3*stEq)]);
    }
  }

  return;
}


void condensedKroneckerProduct(std::vector<mpq_class>& result, const std::vector<mpq_class>& a, const std::vector<mpq_class>& b)
{
  using namespace std;

  vector<mpq_class>::size_type myIndex = 0;

  for (vector<mpq_class>::size_type i = 0; i < a.size(); ++i)
  {
    result[myIndex++] = (a[i]*b[i]);
    for (vector<mpq_class>::size_type j = i+1; j < a.size(); ++j)
    {
      result[myIndex++] = (a[i]*b[j]+a[j]*b[i]);
    }
  }
  return;
}

void sparseMatTimesVec(std::vector<mpq_class>& result, const std::vector<int>& rowInd, const std::vector<int>& colInd, const std::vector<mpq_class>& myCoeffs, const std::vector<mpq_class>& x)
{
  using namespace std;

  std::fill(result.begin(), result.end(), mpq_class(0));

  for(vector<int>::size_type i = 0; i < rowInd.size(); ++i)
  {
    result[rowInd[i]] += myCoeffs[i] * x[colInd[i]];
  }

  return;
}

// END OF PROGRAM
