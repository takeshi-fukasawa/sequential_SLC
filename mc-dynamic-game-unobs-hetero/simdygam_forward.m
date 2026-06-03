% -----------------------------------------------------------------------
%
% SIMDYGAM        Simulates data of state and decision variables
%                 from the steady-state distribution of
%                 a Markov Perfect Equilibrium in a dynamic
%                 game of firms' market entry/exit with
%                 incomplete information
%
% by Victor Aguirregabiria
%
% Converted from Gauss to Matlab by Jason Blevins.
%
% --------------------------------------------------------------------------
%
% FORMAT:
% { aobs, aobs_1, indzexo_obs } = simdygam_forward(nobs,pequil,psteady,mstate)
%
% INPUTS
%     nobs    - Number of simulations (markets)
%
%     pchoice - (nstate x nplayer) matrix with MPE probs of entry
%
%     psteady - (nstate x 1) vector with steady-state distribution of {s[t],a[t-1]}
%
%     mstate  - (nstate x (nplayer+size(z_exoval,2))) matrix with values of state variables {s[t],a[t-1]}.
%              The states in the rows of prob are ordered as follows.
%              Example: z_exo_val=(1|2) and 3 players:
%                       s[t]    a[1,t-1]    a[2,t-1]    a[3,t-1]
%              Row 1:     1           0           0           0
%              Row 2:     1           0           0           1
%              Row 3:     1           0           1           0
%              Row 4:     1           0           1           1
%              Row 5:     1           1           0           0
%              Row 6:     1           1           0           1
%              Row 7:     1           1           1           0
%              Row 8:     1           1           1           1
%              Row 9:     2           0           0           0
%              Row 10:    2           0           0           1
%              Row 11:    2           0           1           0
%              Row 12:    2           0           1           1
%              Row 13:    2           1           0           0
%              Row 14:    2           1           0           1
%              Row 15:    2           1           1           0
%              Row 16:    2           1           1           1
%
% OUTPUTS:
%     aobs   - (nobs x nplayer) matrix with players' choices.
%
%     aobs_1 - (nobs x nplayer) matrix with players' initial states.
%
%     indzexo_obs - (nobs x 1) vector with simulated values of s[t]
%
% ----------------------------------------------------------------------------

function [ aobs, aobs_1, zexo_obs, indobs_vec ] = simdygam_forward(nobs, pchoice, psteady, mstate,ptran,T)

  nplay = size(pchoice, 2);

  % ----------------------------------------------------------------------
  % a. Generating random draws from ergodic distribution of (s[t],a[t-1])
  % ----------------------------------------------------------------------

  indobs_t=random_draw_func(psteady,nobs);%nobs*1

  % --------------------------------------------------------
  % b. Generating random draws for a[t] (given s[t],a[t-1])
  % --------------------------------------------------------
  pchoice = pchoice(indobs_t,:);
  uobs = rand(nobs, nplay);
  aobs = (uobs <= pchoice);

indobs_vec=indobs_t;

%%%%%%%
if T>=2
   indobs_mat=NaN(nobs,T+1);
   indobs_mat(:,1)=indobs_t;

   % ptran: nstate*nstate-dim matrix; Satisfies psteady = ptrana' * psteady
   % ptran(i,j)=Prob(State j|State i)

    for t=2:T+1
       for i=1:nobs
         indobs_t(i)=random_draw_func(reshape(ptran(indobs_t(i),:),[],1),1);
       end
        
       indobs_mat(:,t)=indobs_t;
    end


  indobs_vec=(reshape(indobs_mat(:,1:T),nobs,T))';%T*nobs
  indobs_vec=indobs_vec(:);%(T*nobs)*1

  indobs_vec_plus_1=(reshape(indobs_mat(:,2:T+1),nobs,T))';%T*nobs
  indobs_vec_plus_1=indobs_vec_plus_1(:);

  aobs=mstate(indobs_vec_plus_1,2:nplay+1);

end%if T>=2

  zexo_obs = mstate(indobs_vec, 1);
  aobs_1 = mstate(indobs_vec, 2:nplay+1);

  if size(mstate)>=nplay+2
    zexo_obs=[zexo_obs,mstate(indobs_vec,nplay+2:end)];
  end

end
