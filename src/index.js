'use strict'
/**
 * Created by acastillo on 1/25/16.
 */
function autoassigner(){
    var timeoutTerminated;
    var likelihood;
    var nSolutions=0;
    const DEBUG = false;

    function exploreTreeRec( signals, diaList, indexSignal, indexDia, diaMask, partial){
        //If this happens, we can assign this atom group to this signal
        while(indexDia>=0&&signals[indexSignal]>=diaList[indexDia]){
            //Force a return if the loop time is longer than the given timeout
            const d = new Date();
            if((d.getMilliseconds()-timeStart)>timeout){
                timeoutTerminated=true;
                return;
            }
            //We can speed up it by checking the chemical shift first
            if(diaMask[indexDia]&&isWithinCSRange(indexSignal,indexDia)){
                this.nSteps++;
                const sizePartial = partial[indexSignal].length;
                //Assign the atom group to the signal
                diaMask[indexDia]=false;//Mark this atom group as used
                partial[indexSignal][sizePartial]=indexDia;//Add the atom group index to the assignment list
                signals[indexSignal]-=diaList[indexDia];//Subtract the group from signal integral
                //System.out.println("PartialX "+partial);
                //If this signal is completely assigned, we should verify all the restrictions
                if(signals[indexSignal]==0){
                    let keySum  = accomplishCounts(indexSignal, partial);
                    //System.out.println("Accomplish count: "+keySum);
                    if(keySum!=0){

                        //Verify the restrictions. A good solution should give a high likelihood
                        likelihood = solutionScore(partial, indexSignal, keySum);
                        if(DEBUG) console.log(likelihood+" "+partial);
                        if(DEBUG) console.log("likelihood: "+likelihood);
                        //This is a solution
                        if(likelihood>0){
                            if(DEBUG) console.log(likelihood+" "+partial);
                            //var solution =  {assignment:[],likelihood:likelihood};
                            if(indexSignal==0){//We found a new solution
                                nSolutions++;
                                var solution = {assignment:partial, likelihood:likelihood};
                                //solution.setResult(new JSONArray(partial.toString()), likelihood);
                                //System.out.println(likelihood+" "+partial);
                                //dady.addChild(new SolutionTree(solution));
                                if (solutions.length>=maxSolutions) {
                                    //System.out.println("Solution "+solutions.size()+" "+solutions.last().likelihood);
                                    if (likelihood>solutions[solutions.length-1].likelihood) {
                                        //System.out.println(likelihood+" "+partial);
                                        solutions.pop();
                                        insert(solution, solutions);
                                    }
                                } else {
                                    insert(solution, solutions);
                                }
                                //System.out.println(solutions);
                            }
                            else{
                                //Each new signal that we assign will produce a new level on the tree.
                                indexSignal--;//Lets go forward with the next signal
                                indexDia=diaList.length;
                                while(!diaMask[--indexDia]);
                                //SolutionTree newLevel = new SolutionTree(solution);
                                exploreTreeRec(signals, diaList,indexSignal, indexDia, diaMask, partial);
                                //exploreTreeRec(signals, diaList,indexSignal, indexDia, diaMask, newLevel,partial);
                                //dady.addChild(newLevel);
                                indexSignal++;
                            }
                        }
                    }

                }
                else{
                    //System.out.println("Still assinging "+indexSignal);
                    //It says that the signal should be assigned by combining 2 or more signals
                    const previousIndexDia = indexDia;
                    while(indexDia>0&&!diaMask[--indexDia]);
                    if(indexDia>=0)
                        exploreTreeRec(signals, diaList,indexSignal, indexDia, diaMask, partial);
                    indexDia = previousIndexDia;
                    //return
                }
                //Deallocate this atom group to try the next one.
                indexDia=partial[indexSignal].splice(sizePartial,1);//Add the atom group index to the assignment list
                diaMask[indexDia]=true;//Mark this atom group as available
                //System.out.println("After remove "+partial);
                signals[indexSignal]+=diaList[indexDia];//Subtract the group from signal integral
                likelihoods[indexSignal]=1;//-1
            }//[[4,3,1],[0],[2]]
            indexDia--;
        }
    }

    function insert(element, array) {
        array.splice(locationOf(element, array) + 1, 0, element);
        return array;
    }

    function locationOf(element, array, start, end) {
        start = start || 0;
        end = end || array.length;
        var pivot = parseInt(start + (end - start) / 2, 10);
        if (end-start <= 1 || array[pivot].likelihood === element.likelihood) return pivot;
        if (array[pivot].likelihood < element) {
            return locationOf(element, array, pivot, end);
        } else {
            return locationOf(element, array, start, pivot);
        }
    }

}
