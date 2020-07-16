#include <permon/private/qpimpl.h>
#include <permonpc.h>
#include <permon/private/qppfimpl.h>
#include <permonqpfeti.h>

PetscLogEvent QPT_HomogenizeEq, QPT_EnforceEqByProjector, QPT_EnforceEqByPenalty, QPT_OrthonormalizeEq, QPT_SplitBE;
PetscLogEvent QPT_Dualize, QPT_Dualize_AssembleG, QPT_Dualize_FactorK, QPT_Dualize_PrepareBt, QPT_FetiPrepare, QPT_AllInOne, QPT_RemoveGluingOfDirichletDofs;

static QPPF QPReusedCP = NULL;

/* common tasks during a QP transform - should be called in the beginning of each transform function */
#undef __FUNCT__
#define __FUNCT__ "QPTransformBegin_Private"
static PetscErrorCode QPTransformBegin_Private(PetscErrorCode(*transform)(QP), const char *trname,
    PetscErrorCode(*postSolve)(QP,QP), PetscErrorCode(*postSolveCtxDestroy)(void*),
    QPDuplicateOption opt, QP *qp_inout, QP *child_new, MPI_Comm *comm)
{
  QP child;
  QP qp = *qp_inout;

  PetscFunctionBegin;
  TRY( PetscObjectGetComm((PetscObject)qp,comm) );
  TRY( QPChainGetLast(qp,&qp) );
  TRY( QPSetUpInnerObjects(qp) );
  TRY( QPChainAdd(qp,opt,&child) );
  child->transform = transform;
  TRY( PetscStrcpy(child->transform_name, trname) );
  TRY( QPSetPC(child,qp->pc) );
  if (qp->changeListener) TRY( (*qp->changeListener)(qp) );

  child->postSolve = QPDefaultPostSolve;
  if (postSolve) child->postSolve = postSolve;
  child->postSolveCtxDestroy = postSolveCtxDestroy;

  TRY( FllopPetscObjectInheritPrefix((PetscObject)child,(PetscObject)qp,NULL) );

  *qp_inout = qp;
  *child_new = child;
  TRY( PetscInfo5(qp,"QP %x (#%d) transformed by %s to QP %x (#%d)\n",qp,qp->id,trname,child,child->id) );
  PetscFunctionReturn(0);
}

#define QPTransformBegin(transform,postSolve,postSolveCtxDestroy,opt,qp,child,comm) QPTransformBegin_Private((PetscErrorCode(*)(QP))transform,__FUNCT__,\
    (PetscErrorCode(*)(QP,QP))postSolve, (PetscErrorCode(*)(void*))postSolveCtxDestroy,\
    opt,qp,child,comm)

#undef __FUNCT__
#define __FUNCT__ "QPDefaultPostSolve"
PetscErrorCode QPDefaultPostSolve(QP child,QP parent)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(child,QP_CLASSID,1);
  PetscValidHeaderSpecific(parent,QP_CLASSID,2);
  TRY( VecCopy(child->x,parent->x) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTEnforceEqByProjectorPostSolve_Private"
static PetscErrorCode QPTEnforceEqByProjectorPostSolve_Private(QP child,QP parent)
{
  Mat A = parent->A;
  Vec b = parent->b;
  Vec x = parent->x;
  Vec lambda_E  = parent->lambda_E;
  Vec Bt_lambda = parent->Bt_lambda;
  Vec r = parent->xwork;
  PetscBool skip_lambda_E=PETSC_TRUE, skip_Bt_lambda=PETSC_TRUE;
  PetscBool inherit_eq_multipliers=PETSC_TRUE;

  PetscFunctionBegin;
  TRY( QPDefaultPostSolve(child,parent) );

  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_project_inherit_eq_multipliers",&inherit_eq_multipliers,NULL) );

  if (child->BE && inherit_eq_multipliers) {
    TRY( VecIsInvalidated(child->lambda_E,&skip_lambda_E) );
    TRY( VecIsInvalidated(child->Bt_lambda,&skip_Bt_lambda) );
  }

  TRY( MatMult(A, x, r) );                                                      /* r = A*x */
  TRY( VecAYPX(r,-1.0,b) );                                                     /* r = b - r */

  if (!skip_lambda_E) {
    /* lambda_E1 = lambda_E2 + (B*B')\B*(b-A*x) */
    TRY( QPPFApplyHalfQ(parent->pf,r,lambda_E) );
    TRY( VecAXPY(lambda_E,1.0,child->lambda_E) );
  }
  if (!skip_Bt_lambda) {
    /* (B'*lambda)_1 = (B'*lambda)_2 + B'*(B*B')\B*(b-A*x) */
    TRY( QPPFApplyQ(parent->pf,r,Bt_lambda)  );
    TRY( VecAXPY(Bt_lambda,1.0,child->Bt_lambda) );
  }
  PetscFunctionReturn(0);
}

typedef struct {
  Mat  P;
  PC   pc;
  Vec  work;
  PetscBool symmetric;
} PC_QPTEnforceEqByProjector;

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_QPTEnforceEqByProjector"
static PetscErrorCode PCDestroy_QPTEnforceEqByProjector(PC pc)
{
  PC_QPTEnforceEqByProjector* ctx;

  PetscFunctionBegin;
  TRY( PCShellGetContext(pc,(void**)&ctx) );
  TRY( MatDestroy(&ctx->P) );
  TRY( PCDestroy(&ctx->pc) );
  TRY( VecDestroy(&ctx->work) );
  TRY( PetscFree(ctx) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_QPTEnforceEqByProjector"
static PetscErrorCode PCApply_QPTEnforceEqByProjector(PC pc,Vec x,Vec y)
{
  PC_QPTEnforceEqByProjector* ctx;

  PetscFunctionBegin;
  TRY( PCShellGetContext(pc,(void**)&ctx) );
  TRY( PCApply(ctx->pc,x,ctx->work) );
  TRY( MatMult(ctx->P,ctx->work,y) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_QPTEnforceEqByProjector_Symmetric"
static PetscErrorCode PCApply_QPTEnforceEqByProjector_Symmetric(PC pc,Vec x,Vec y)
{
  PC_QPTEnforceEqByProjector* ctx;

  PetscFunctionBegin;
  TRY( PCShellGetContext(pc,(void**)&ctx) );
  TRY( MatMult(ctx->P,x,y) );
  TRY( PCApply(ctx->pc,y,ctx->work) );
  TRY( MatMult(ctx->P,ctx->work,y) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_QPTEnforceEqByProjector_None"
static PetscErrorCode PCApply_QPTEnforceEqByProjector_None(PC pc,Vec x,Vec y)
{
  PetscFunctionBegin;
  TRY( VecCopy(x,y) );
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PCSetUp_QPTEnforceEqByProjector"
static PetscErrorCode PCSetUp_QPTEnforceEqByProjector(PC pc)
{
  PC_QPTEnforceEqByProjector* ctx;
  PetscBool flg;
  PCDualType type;
  PetscBool none;

  PetscFunctionBegin;
  TRY( PCShellGetContext(pc,(void**)&ctx) );

  none = PETSC_FALSE;
  TRY( PetscObjectTypeCompare((PetscObject)ctx->pc,PCNONE,&none) );
  if (!none) {
    TRY( PetscObjectTypeCompare((PetscObject)ctx->pc,PCDUAL,&flg) );
    if (flg) {
      TRY( PCDualGetType(ctx->pc,&type) );
      if (type == PC_DUAL_NONE) none = PETSC_TRUE;
    }
  }

  if (none) {
    TRY( PCShellSetApply(pc,PCApply_QPTEnforceEqByProjector_None) );
  } else if (ctx->symmetric) {
    TRY( PCShellSetApply(pc,PCApply_QPTEnforceEqByProjector_Symmetric) );
  } else {
    TRY( PCShellSetApply(pc,PCApply_QPTEnforceEqByProjector) );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCCreate_QPTEnforceEqByProjector"
PetscErrorCode PCCreate_QPTEnforceEqByProjector(PC pc_orig,Mat P,PetscBool symmetric,PC *pc_new)
{
  PC_QPTEnforceEqByProjector *ctx;
  PC pc;

  PetscFunctionBegin;
  TRY( PetscNew(&ctx) );
  ctx->symmetric = symmetric;
  ctx->P = P;
  ctx->pc = pc_orig;
  TRY( PetscObjectReference((PetscObject)P) );
  TRY( PetscObjectReference((PetscObject)pc_orig) );
  TRY( MatCreateVecs(P,&ctx->work,NULL) );

  TRY( PCCreate(PetscObjectComm((PetscObject)pc_orig),&pc) );
  TRY( PCSetType(pc,PCSHELL) );
  TRY( PCShellSetName(pc,"QPTEnforceEqByProjector") );
  TRY( PCShellSetContext(pc,ctx) );
  TRY( PCShellSetApply(pc,PCApply_QPTEnforceEqByProjector) );
  TRY( PCShellSetDestroy(pc,PCDestroy_QPTEnforceEqByProjector) );
  TRY( PCShellSetSetUp(pc,PCSetUp_QPTEnforceEqByProjector) );

  *pc_new = pc;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTEnforceEqByProjector"
PetscErrorCode QPTEnforceEqByProjector(QP qp)
{
  MPI_Comm         comm;
  Mat              P=NULL;
  Mat              newA=NULL, A_arr[3];
  Vec              newb=NULL;
  QP               child;
  PetscBool        eqonly=PETSC_FALSE;
  PetscBool        pc_symmetric=PETSC_FALSE;
  PetscBool        inherit_box_multipliers=PETSC_FALSE;

  FllopTracedFunctionBegin;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);

  TRY( QPChainGetLast(qp,&qp) );
  if (!qp->BE) {
    TRY( PetscInfo(qp, "no lin. eq. con. matrix specified ==> nothing to enforce\n") );
    TRY( VecDestroy(&qp->cE) );
    PetscFunctionReturn(0);
  }

  FllopTraceBegin;
  if (qp->cE) {
    TRY( PetscInfo(qp, "nonzero lin. eq. con. RHS prescribed ==> automatically calling QPTHomogenizeEq\n") );
    TRY( QPTHomogenizeEq(qp) );
    TRY( QPChainGetLast(qp,&qp) );
  }

  TRY( PetscLogEventBegin(QPT_EnforceEqByProjector,qp,0,0,0) );
  TRY( QPTransformBegin(QPTEnforceEqByProjector, QPTEnforceEqByProjectorPostSolve_Private,NULL, QP_DUPLICATE_DO_NOT_COPY,&qp,&child,&comm) );
  TRY( QPAppendOptionsPrefix(child,"proj_") );

  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_project_pc_symmetric",&pc_symmetric,NULL) );
  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_project_inherit_box_multipliers",&inherit_box_multipliers,NULL) );

  eqonly = !(qp->BI || qp->qpc);
  if (eqonly) {
    TRY( PetscInfo(qp, "only lin. eq. con. were prescribed ==> they are now eliminated\n") );
    TRY( QPSetEq(  child, NULL, NULL) );
  } else {
    TRY( PetscInfo(qp, "NOT only lin. eq. con. prescribed ==> lin. eq. con. are NOT eliminated\n") );
    TRY( QPSetQPPF(child, qp->pf) );
    TRY( QPSetEq(  child, qp->BE, qp->cE) );
  }
  TRY( QPSetIneq(child, qp->BI, qp->cI) );
  TRY( QPSetRhs( child, qp->b) );

  if (inherit_box_multipliers) {
    TRY( QPSetQPC(child, qp->qpc) );
  } else {
    /* TODO: generalize for other QPC */
    Vec lb,ub;
    IS is;
    TRY( QPGetBox(qp,&is,&lb,&ub) );
    TRY( QPSetBox(child,is,lb,ub) );
  }

  TRY( QPPFCreateP(qp->pf,&P) );
  if (eqonly) {
    /* newA = P*A */
    A_arr[0]=qp->A; A_arr[1]=P;
    TRY( MatCreateProd(comm,2,A_arr,&newA) );
  } else {
    /* newA = P*A*P */
    A_arr[0]=P; A_arr[1]=qp->A; A_arr[2]=P;
    TRY( MatCreateProd(comm,3,A_arr,&newA) );
  }
  TRY( QPSetOperator(child,newA) );
  TRY( QPSetOperatorNullSpace(child,qp->R) );
  TRY( MatDestroy(&newA) );

  /* newb = P*b */
  TRY( VecDuplicate(qp->b, &newb) );
  TRY( MatMult(P, qp->b, newb) );
  if (FllopDebugEnabled) {
    PetscReal norm1, norm2;
    TRY( VecNorm(qp->b, NORM_2, &norm1) );
    TRY( VecNorm(newb, NORM_2, &norm2) );
    TRY( FllopDebug2("\n    ||b||\t= %.12e  b = b_bar\n    ||Pb||\t= %.12e\n",norm1,norm2) );
  }
  TRY( QPSetRhs(child, newb) );
  TRY( VecDestroy(&newb) );

  /* create special preconditioner pc_child = P * pc_parent */
  {
    PC pc_parent,pc_child;

    TRY( QPGetPC(qp,&pc_parent) );
    TRY( PCCreate_QPTEnforceEqByProjector(pc_parent,P,pc_symmetric,&pc_child) );
    TRY( QPSetPC(child,pc_child) );
    TRY( PCDestroy(&pc_child) );
  }

  TRY( QPSetWorkVector(child,qp->xwork) );

  TRY( MatDestroy(&P) );
  TRY( PetscLogEventEnd(QPT_EnforceEqByProjector,qp,0,0,0) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTEnforceEqByPenalty_PostSolve_Private"
static PetscErrorCode QPTEnforceEqByPenalty_PostSolve_Private(QP child,QP parent)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

//TODO allow to set multiple of max eigenvalue
#undef __FUNCT__
#define __FUNCT__ "QPTEnforceEqByPenalty"
PetscErrorCode QPTEnforceEqByPenalty(QP qp, PetscReal rho_user, PetscBool rho_direct)
{
  MPI_Comm         comm;
  Mat              newA=NULL;
  Vec              newb=NULL;
  QP               child;
  PetscReal        rho,maxeig;
  PetscReal        maxeig_tol=PETSC_DECIDE;
  PetscInt         maxeig_iter=PETSC_DECIDE;

  FllopTracedFunctionBegin;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  PetscValidLogicalCollectiveReal(qp,rho_user,2);
  PetscValidLogicalCollectiveBool(qp,rho_direct,3);

  if (!rho_user) {
    TRY( PetscInfo(qp, "penalty=0.0 ==> no effect, returning...\n") );
    PetscFunctionReturn(0);
  }

  TRY( QPChainGetLast(qp,&qp) );

  if (!qp->BE) {
    TRY( PetscInfo(qp, "no lin. eq. con. matrix specified ==> nothing to enforce\n") );
    TRY( VecDestroy(&qp->cE) );
    PetscFunctionReturn(0);
  }

  FllopTraceBegin;
  if (qp->cE) {
    PetscBool flg=PETSC_FALSE;
    TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_homogenize_eq_always",&flg,NULL) );
    if (flg) {
      TRY( PetscInfo(qp, "nonzero lin. eq. con. RHS prescribed and -qpt_homogenize_eq_always set to true ==> automatically calling QPTHomogenizeEq\n") );
      TRY( QPTHomogenizeEq(qp) );
      TRY( QPChainGetLast(qp,&qp) );
    }
  }

  TRY( PetscOptionsGetReal(NULL,NULL,"-qpt_penalize_maxeig_tol",&maxeig_tol,NULL) );
  TRY( PetscOptionsGetInt(NULL,NULL,"-qpt_penalize_maxeig_iter",&maxeig_iter,NULL) );

  if (!rho_direct) {
    TRY( MatGetMaxEigenvalue(qp->A, NULL, &maxeig, maxeig_tol, maxeig_iter) );
    rho = rho_user * maxeig;
  } else {
    rho = rho_user;
  }

  if (rho < 0) {
    FLLOP_SETERRQ(PetscObjectComm((PetscObject)qp),PETSC_ERR_ARG_WRONG,"rho must be nonnegative");
  }
  TRY( PetscInfo1(qp, "using penalty = real %.12e\n", rho) );

  TRY( PetscLogEventBegin(QPT_EnforceEqByPenalty,qp,0,0,0) );
  TRY( QPTransformBegin(QPTEnforceEqByPenalty, QPTEnforceEqByPenalty_PostSolve_Private,NULL, QP_DUPLICATE_DO_NOT_COPY,&qp,&child,&comm) );
  TRY( QPAppendOptionsPrefix(child,"pnlt_") );

  TRY( QPSetEq(  child, NULL, NULL) );
  TRY( QPSetQPC(child, qp->qpc) );
  TRY( QPSetIneq(child, qp->BI, qp->cI) );
  TRY( QPSetRhs(child, qp->b) );
  TRY( QPSetInitialVector(child, qp->x) );

  /* newA = A + rho*BE'*BE */
  TRY( MatCreatePenalized(qp,rho,&newA) );
  TRY( QPSetOperator(child,newA) );
  TRY( QPSetOperatorNullSpace(child,qp->R) );
  TRY( MatDestroy(&newA) );

  /* newb = b + rho*BE'*c */
  if (qp->c) {
    TRY( VecDuplicate(qp->b,&newb) );
    TRY( MatMultTranspose(qp->BE,qp->c,newb) );
    TRY( VecAYPX(newb,rho,qp->b) );
    TRY( QPSetRhs(child,newb) );
  }
  
  TRY( QPSetIneqMultiplier(       child,qp->lambda_I) );
  TRY( QPSetWorkVector(child,qp->xwork) );

  TRY( PetscLogEventEnd(QPT_EnforceEqByPenalty,qp,0,0,0) );
  PetscFunctionReturnI(0);
}


#undef __FUNCT__
#define __FUNCT__ "QPTHomogenizeEqPostSolve_Private"
static PetscErrorCode QPTHomogenizeEqPostSolve_Private(QP child,QP parent)
{
  Vec xtilde = (Vec) child->postSolveCtx;

  PetscFunctionBegin;
  /* x_parent = x_child + xtilde */
  TRY( VecWAXPY(parent->x,1.0,child->x,xtilde) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTHomogenizeEqPostSolveCtxDestroy_Private"
static PetscErrorCode QPTHomogenizeEqPostSolveCtxDestroy_Private(void *ctx)
{
  Vec xtilde = (Vec) ctx;

  PetscFunctionBegin;
  TRY( VecDestroy(&xtilde) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTHomogenizeEq"
PetscErrorCode QPTHomogenizeEq(QP qp)
{
  MPI_Comm          comm;
  QP                child;
  Vec               b_bar, cineq, lb, ub, lbnew, ubnew, xtilde, xtilde_sub;
  IS                is;

  FllopTracedFunctionBegin;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  lb = NULL; lbnew=NULL;
  ub = NULL; ubnew=NULL;
  is = NULL;

  TRY( QPChainGetLast(qp,&qp) );
  if (!qp->cE) {
    TRY( PetscInfo(qp, "lin. eq. con. already homogenous ==> nothing to homogenize\n") );
    PetscFunctionReturn(0);
  }

  FllopTraceBegin;
  TRY( PetscLogEventBegin(QPT_HomogenizeEq,qp,0,0,0) );
  TRY( QPTransformBegin(QPTHomogenizeEq, QPTHomogenizeEqPostSolve_Private,QPTHomogenizeEqPostSolveCtxDestroy_Private, QP_DUPLICATE_COPY_POINTERS,&qp,&child,&comm) );

  /* A, R remain the same */

  TRY( VecDuplicate(qp->x,&xtilde) );
  TRY( QPPFApplyHalfQTranspose(qp->pf,qp->cE,xtilde) );                         /* xtilde = BE'*inv(BE*BE')*cE */

  TRY( VecDuplicate(qp->b, &b_bar) );
  TRY( MatMult(qp->A, xtilde, b_bar) );
  TRY( VecAYPX(b_bar, -1.0, qp->b) );                                           /* b_bar = b - A*xtilde */
  TRY( QPSetRhs(child, b_bar) );
  TRY( VecDestroy(&b_bar) );

  if (FllopDebugEnabled) {
    PetscReal norm1,norm2,norm3,norm4;
    TRY( VecNorm(qp->cE, NORM_2, &norm1) );
    TRY( VecNorm(xtilde, NORM_2, &norm2) );
    TRY( VecNorm(qp->b,  NORM_2, &norm3) );
    TRY( VecNorm(child->b,NORM_2, &norm4) );
    TRY( FllopDebug4("\n"
        "    ||ceq||\t= %.12e  ceq = e\n"
        "    ||xtilde||\t= %.12e  xtilde = Beq'*inv(Beq*Beq')*ceq\n"
        "    ||b||\t= %.12e  b = d\n"
        "    ||b_bar||\t= %.12e  b_bar = b-A*xtilde\n",norm1,norm2,norm3,norm4) );
  }

  TRY( QPSetQPPF(child, qp->pf) );
  TRY( QPSetEq(child,qp->BE,NULL) );                                            /* cE is eliminated */

  cineq = NULL;
  if (qp->cI) {
    TRY( VecDuplicate(qp->cI, &cineq) );
    TRY( MatMult(qp->BI, xtilde, cineq) );
    TRY( VecAYPX(cineq, -1.0, qp->cI) );                                        /* cI = cI - BI*xtilde */
  }
  TRY( QPSetIneq(child, qp->BI, cineq) );
  TRY( VecDestroy(&cineq) );

  TRY( QPGetBox(qp, &is, &lb, &ub) );
  if (is) {
    TRY( VecGetSubVector(xtilde, is, &xtilde_sub) );
  } else {
    xtilde_sub = xtilde;
  }

  if (lb) {
    TRY( VecDuplicate(lb, &lbnew) );
    TRY( VecWAXPY(lbnew, -1.0, xtilde_sub, lb) );                                /* lb = lb - xtilde */
  }

  if (ub) {
    TRY( VecDuplicate(ub, &ubnew) );
    TRY( VecWAXPY(ubnew, -1.0, xtilde_sub, ub) );                                /* ub = ub - xtilde */
  }

  if (is) {
    TRY( VecRestoreSubVector(xtilde, is, &xtilde_sub) );
  }

  TRY( QPSetBox(child, is, lbnew, ubnew) );
  TRY( VecDestroy(&lbnew) );
  TRY( VecDestroy(&ubnew) );

  TRY( VecDestroy(&child->x) );

  child->postSolveCtx = xtilde;
  TRY( PetscLogEventEnd(QPT_HomogenizeEq,qp,0,0,0) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolve_QPTOrthonormalizeEq"
static PetscErrorCode QPTPostSolve_QPTOrthonormalizeEq(QP child,QP parent)
{
  Mat T = (Mat) child->postSolveCtx;
  PetscBool skip_lambda_E=PETSC_TRUE,skip_Bt_lambda=PETSC_TRUE;

  PetscFunctionBegin;
  /* B_p'*lambda_p = B_c'*lambda_c is done implicitly */

  //TRY( VecIsInvalidated(child->lambda_E,&skip_lambda_E) );
  //TRY( VecIsInvalidated(child->Bt_lambda,&skip_Bt_lambda) );

  if (!skip_lambda_E && T->ops->multtranspose) {
    /* lambda_E_p = T'*lambda_E_c */
    TRY( MatMultTranspose(T,child->lambda_E,parent->lambda_E) );
  } else if (!skip_Bt_lambda) {
    //TODO this seems to be inaccurate for explicit orthonormalization
    /* lambda_E_p = (G*G')\G(G'*lambda_c) where G'*lambda_c = child->Bt_lambda */
    TRY( QPPFApplyHalfQ(child->pf,child->Bt_lambda,parent->lambda_E) );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolveDestroy_QPTOrthonormalizeEq"
static PetscErrorCode QPTPostSolveDestroy_QPTOrthonormalizeEq(void *ctx)
{
  Mat T = (Mat) ctx;

  PetscFunctionBegin;
  TRY( MatDestroy(&T) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTOrthonormalizeEq"
PetscErrorCode QPTOrthonormalizeEq(QP qp,MatOrthType type,MatOrthForm form)
{
  MPI_Comm          comm;
  QP                child;
  Mat               BE,TBE,T;
  Vec               cE,TcE;

  FllopTracedFunctionBegin;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);

  TRY( QPChainGetLast(qp,&qp) );
  if (!qp->BE) {
    TRY( PetscInfo(qp, "no lin. eq. con. ==> nothing to orthonormalize\n") );
    PetscFunctionReturn(0);
  }
  if (type == MAT_ORTH_NONE) {
    TRY( PetscInfo(qp, "MAT_ORTH_NONE ==> skipping orthonormalization\n") );
    PetscFunctionReturn(0);
  }

  FllopTraceBegin;
  if (type == MAT_ORTH_INEXACT && qp->cE) {
    TRY( PetscInfo(qp, "MAT_ORTH_INEXACT and nonzero lin. eq. con. RHS prescribed ==> automatically calling QPTHomogenizeEq\n") );
    TRY( QPTHomogenizeEq(qp) );
    TRY( QPChainGetLast(qp,&qp) );
  }

  TRY( PetscLogEventBegin(QPT_OrthonormalizeEq,qp,0,0,0) );
  TRY( QPTransformBegin(QPTOrthonormalizeEq, QPTPostSolve_QPTOrthonormalizeEq,QPTPostSolveDestroy_QPTOrthonormalizeEq, QP_DUPLICATE_COPY_POINTERS,&qp,&child,&comm) );

  TRY( QPGetEq(qp, &BE, &cE) );
  TRY( MatOrthRows(BE, type, form, &TBE, &T) );

  {
    const char *name;
    TRY( PetscObjectGetName((PetscObject)BE,&name) );
    TRY( PetscObjectSetName((PetscObject)TBE,name) );
    TRY( MatPrintInfo(T) );
    TRY( MatPrintInfo(TBE) );
  }

  TcE=NULL;
  if (cE) {
    if (type == MAT_ORTH_IMPLICIT || type == MAT_ORTH_INEXACT) {
      TRY( PetscObjectReference((PetscObject)cE) );
      TcE = cE;
    } else {
      TRY( VecDuplicate(cE,&TcE) );
      TRY( MatMult(T,cE,TcE) );
    }
  }

  TRY( QPSetQPPF(child,NULL) );
  TRY( QPSetEq(child,TBE,TcE) );                /* QPPF re-created in QPSetEq */
  TRY( QPSetEqMultiplier(child,NULL) );         /* lambda_E will be re-created in QPSetUp */

  if (type == MAT_ORTH_IMPLICIT) {
    /* inherit GGtinv to avoid extra GGt factorization - dirty way */
    TRY( QPPFSetUp(qp->pf) );
    child->pf->GGtinv = qp->pf->GGtinv;
    child->pf->Gt = qp->pf->Gt;
    TRY( PetscObjectReference((PetscObject)qp->pf->GGtinv) );
    TRY( PetscObjectReference((PetscObject)qp->pf->Gt) );
    TRY( QPPFSetUp(child->pf) );
  } else if (type == MAT_ORTH_INEXACT) {
    TRY( PetscObjectSetOptionsPrefix((PetscObject)child->pf,"inexact_") );
    TRY( PetscObjectCompose((PetscObject)child->pf,"exact",(PetscObject)qp->pf) );
  }

  TRY( MatDestroy(&TBE) );
  TRY( VecDestroy(&TcE) );
  child->postSolveCtx = T;
  TRY( PetscLogEventEnd(QPT_OrthonormalizeEq,qp,0,0,0) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTOrthonormalizeEqFromOptions"
PetscErrorCode QPTOrthonormalizeEqFromOptions(QP qp)
{
  MatOrthType eq_orth_type=MAT_ORTH_NONE;
  MatOrthForm eq_orth_form=MAT_ORTH_FORM_IMPLICIT;
  QP last;

  PetscFunctionBeginI;
  TRY( QPChainGetLast(qp,&last) );
  _fllop_ierr = PetscObjectOptionsBegin((PetscObject)last);CHKERRQ(_fllop_ierr);
  TRY( PetscOptionsEnum("-qp_E_orth_type","type of eq. con. orthonormalization","QPTOrthonormalizeEq",MatOrthTypes,(PetscEnum)eq_orth_type,(PetscEnum*)&eq_orth_type,NULL) );
  TRY( PetscOptionsEnum("-qp_E_orth_form","form of eq. con. orthonormalization","QPTOrthonormalizeEq",MatOrthForms,(PetscEnum)eq_orth_form,(PetscEnum*)&eq_orth_form,NULL) );
  TRY( PetscInfo1(qp, "-qp_E_orth_type %s\n",MatOrthTypes[eq_orth_type]) );
  TRY( PetscInfo1(qp, "-qp_E_orth_form %s\n",MatOrthForms[eq_orth_form]) );
  TRY( QPTOrthonormalizeEq(last,eq_orth_type,eq_orth_form) );
  _fllop_ierr = PetscOptionsEnd();CHKERRQ(_fllop_ierr);
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTDualizeViewBSpectra_Private"
static PetscErrorCode QPTDualizeViewBSpectra_Private(Mat B)
{
  PetscFunctionBegin;
  PetscInt j,Mn;
  Mat Bt,BBt,Bn;
  PetscReal norm;
  PetscBool flg;
  const char *name;

  TRY( PetscObjectTypeCompareAny((PetscObject)B,&flg,MATNEST,MATNESTPERMON,"") );
  if (flg) {
    TRY( MatNestGetSize(B,&Mn,NULL) );
    for (j=0; j<Mn; j++) {
      TRY( MatNestGetSubMat(B,j,0,&Bn) );
      TRY( QPTDualizeViewBSpectra_Private(Bn) );
    }
  }
  TRY( PetscObjectGetName((PetscObject)B,&name) );
  TRY( MatCreateTranspose(B,&Bt) );
  TRY( MatCreateNormal(Bt,&BBt) );
  TRY( MatGetMaxEigenvalue(BBt,NULL,&norm,1e-6,10) );
  TRY( PetscPrintf(PetscObjectComm((PetscObject)B),"||%s * %s'|| = %.12e (10 power method iterations)\n",name,name,norm) );
  TRY( MatDestroy(&BBt) );
  TRY( MatDestroy(&Bt) );
  PetscFunctionReturn(0);
}

#define QPTDualizeView_Private_SetName(mat,matname) if (mat && !((PetscObject)mat)->name) TRY( PetscObjectSetName((PetscObject)mat,matname) )

#undef __FUNCT__
#define __FUNCT__ "QPTDualizeView_Private"
static PetscErrorCode QPTDualizeView_Private(QP qp, QP child)
{
  MPI_Comm comm;
  Mat F,Kplus,K,Kreg,B,Bt,R,G;
  Vec c,d,e,f,lb,lambda,u;

  PetscFunctionBegin;
  TRY( PetscObjectGetComm((PetscObject)qp, &comm) );

  F = child->A;
  TRY( PetscObjectQuery((PetscObject)F,"Kplus",(PetscObject*)&Kplus) );
  TRY( PetscObjectQuery((PetscObject)F,"B",(PetscObject*)&B) );
  TRY( PetscObjectQuery((PetscObject)F,"Bt",(PetscObject*)&Bt) );
  
  K = qp->A;
  R = qp->R;
  G = child->BE;
  c = qp->c;
  d = child->b;
  e = child->cE;
  f = qp->b;
  TRY( QPGetBox(child,NULL,&lb,NULL) );
  lambda = child->x;
  u = qp->x;

  if (Kplus) {
    Mat Kplus_inner;
    TRY( PetscObjectQuery((PetscObject)Kplus,"Kplus",(PetscObject*)&Kplus_inner) );
    if (Kplus_inner) Kplus = Kplus_inner;
    TRY( MatInvGetRegularizedMat(Kplus,&Kreg) );
  }

  QPTDualizeView_Private_SetName(K,     "K");
  QPTDualizeView_Private_SetName(Kreg,  "Kreg");
  QPTDualizeView_Private_SetName(Kplus, "Kplus");
  QPTDualizeView_Private_SetName(R,     "R");
  QPTDualizeView_Private_SetName(u,     "u");
  QPTDualizeView_Private_SetName(f,     "f");
  QPTDualizeView_Private_SetName(B,     "B");
  QPTDualizeView_Private_SetName(Bt,    "Bt");
  QPTDualizeView_Private_SetName(c,     "c");
  QPTDualizeView_Private_SetName(F,     "F");
  QPTDualizeView_Private_SetName(lambda,"lambda");
  QPTDualizeView_Private_SetName(d,     "d");
  QPTDualizeView_Private_SetName(G,     "G");
  QPTDualizeView_Private_SetName(e,     "e");
  QPTDualizeView_Private_SetName(lb,    "lb");

  if (FllopObjectInfoEnabled && !PetscPreLoadingOn) {
    TRY( PetscPrintf(comm, "*** "__FUNCT__":\n") );
    if (K)      TRY( MatPrintInfo(K) );
    if (Kreg)   TRY( MatPrintInfo(Kreg) );
    if (Kplus)  TRY( MatPrintInfo(Kplus) );
    if (R)      TRY( MatPrintInfo(R) );
    if (u)      TRY( VecPrintInfo(u) );
    if (f)      TRY( VecPrintInfo(f) );
    if (B)      TRY( MatPrintInfo(B) );
    if (Bt)     TRY( MatPrintInfo(Bt) );
    if (c)      TRY( VecPrintInfo(c) );
    if (F)      TRY( MatPrintInfo(F) );
    if (lambda) TRY( VecPrintInfo(lambda) );
    if (d)      TRY( VecPrintInfo(d) );
    if (G)      TRY( MatPrintInfo(G) );
    if (e)      TRY( VecPrintInfo(e) );
    if (lb)     TRY( VecPrintInfo(lb) );

    TRY( PetscPrintf(comm, "***\n\n") );
  }

  if (FllopDebugEnabled) {
    PetscReal norm1=0.0, norm2=0.0, norm3=0.0, norm4=0.0;
    if (f) TRY( VecNorm(f, NORM_2, &norm1) );
    if (c) TRY( VecNorm(c, NORM_2, &norm2) );
    if (d) TRY( VecNorm(d, NORM_2, &norm3) );
    if (e) TRY( VecNorm(e, NORM_2, &norm4) );
    TRY( FllopDebug4("\n"
        "    ||f||\t= %.12e\n"
        "    ||c||\t= %.12e\n"
        "    ||d||\t= %.12e  d = B*Kplus*f-c \n"
        "    ||e||\t= %.12e  e = R'*f \n",norm1,norm2,norm3,norm4) );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTDualizePostSolve_Private"
static PetscErrorCode QPTDualizePostSolve_Private(QP child,QP parent)
{
    Mat F          = child->A;
    Vec alpha      = child->lambda_E;
    Mat Kplus;
    Vec f          = parent->b;
    Vec u          = parent->x;
    Vec tprim      = parent->xwork;
    PetscBool flg;

    PetscFunctionBegin;
    TRY( PetscObjectQuery((PetscObject)F,"Kplus",(PetscObject*)&Kplus) );
    FLLOP_ASSERT(Kplus,"Kplus != NULL");

    /* copy lambda back to lambda_E and lambda_I */
    if (parent->BE && parent->BI) {
      PetscInt i;
      IS rows[2];
      Vec as[2]={parent->lambda_E,parent->lambda_I};
      Vec a[2];

      TRY( MatNestGetISs(parent->B,rows,NULL) );
      for (i=0; i<2; i++) {
        if (as[i]) {
          TRY( VecGetSubVector(parent->lambda,rows[i],&a[i]) );
          TRY( VecCopy(a[i],as[i]) );
          TRY( VecRestoreSubVector(parent->lambda,rows[i],&a[i]) );
        }
      }
    }

    /* u = Kplus*(f-B'*lambda) */
    TRY( MatMultTranspose(parent->B, parent->lambda, tprim) );
    TRY( VecAYPX(tprim, -1.0, f) );
    TRY( MatMult(Kplus, tprim, u) );

    TRY( VecIsInvalidated(alpha,&flg) );
    if (flg) {
      /* compute alpha = (G*G')\G(G'*alpha) where G'*alpha = child->Bt_lambda */
      TRY( VecIsInvalidated(child->Bt_lambda,&flg) );
      if (flg) {
        TRY( MatMultTranspose(child->B,child->lambda,child->Bt_lambda) );
      }
      TRY( QPPFApplyHalfQ(child->pf,child->Bt_lambda,alpha) );
    }

    /* u = u - R*alpha */
    TRY( MatMult(parent->R, alpha, tprim) );
    TRY( VecAXPY(u, -1.0, tprim) );
    PetscFunctionReturn(0);
}

//TODO this a prototype, integrate to API
#undef __FUNCT__
#define __FUNCT__ "MatTransposeMatMult_R_Bt"
static PetscErrorCode MatTransposeMatMult_R_Bt(Mat R, Mat Bt, Mat *G_new)
{
  Mat       G,Gt;
  PetscBool flg;
  PetscBool G_explicit = PETSC_TRUE;

  PetscFunctionBeginI;
  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_dualize_G_explicit",&G_explicit,NULL) );

  if (!G_explicit) {
    //NOTE doesn't work with redundant GGt_inv & redundancy of R is problematic due to its structure
    Mat G_arr[2];
    Mat Rt;

    TRY( PermonMatTranspose(R,MAT_TRANSPOSE_CHEAPEST,&Rt) );

    TRY( MatCreateTimer(Rt,&G_arr[1]) );
    TRY( MatCreateTimer(Bt,&G_arr[0]) );
    
    TRY( MatCreateProd(PetscObjectComm((PetscObject)Bt), 2, G_arr, &G) );
    TRY( MatCreateTimer(G,&G) );
    
    TRY( PetscObjectCompose((PetscObject)G,"Bt",(PetscObject)Bt) );
    TRY( PetscObjectCompose((PetscObject)G,"R",(PetscObject)R) );
    TRY( PetscObjectSetName((PetscObject)G,"G") );

    TRY( MatDestroy(&Rt) );
    *G_new = G;
    PetscFunctionReturn(0);
  }

  TRY( PetscObjectTypeCompareAny((PetscObject)Bt,&flg,MATNEST,MATNESTPERMON,"") );
  if (!flg) {
    TRY( PetscObjectTypeCompare((PetscObject)Bt,MATEXTENSION,&flg) );
    //TRY( MatTransposeMatMultWorks(R,Bt,&flg) );
    if (flg) {
      TRY( MatTransposeMatMult(R,Bt,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&G) );
    } else {
      TRY( PetscPrintf(PetscObjectComm((PetscObject)Bt), "WARNING: MatTransposeMatMult not applicable, falling back to MatMatMultByColumns\n") );
      TRY( MatTransposeMatMultByColumns(Bt,R,PETSC_TRUE,&Gt) );
      TRY( PermonMatTranspose(Gt,MAT_TRANSPOSE_CHEAPEST,&G) );
      TRY( MatDestroy(&Gt) );
    }
  } else {
    //TODO make this MATNESTPERMON's method
    PetscInt Mn,Nn,J;
    Mat BtJ;
    Mat *mats_G;
    TRY( MatNestGetSize(Bt,&Mn,&Nn) );
    FLLOP_ASSERT(Mn==1,"Mn==1");
    TRY( PetscMalloc(Nn*sizeof(Mat),&mats_G) );
    for (J=0; J<Nn; J++) {
      TRY( MatNestGetSubMat(Bt,0,J,&BtJ) );
      TRY( MatTransposeMatMult_R_Bt(R,BtJ,&mats_G[J]) );
    }
    TRY( MatCreateNestPermon(PetscObjectComm((PetscObject)Bt),1,NULL,Nn,NULL,mats_G,&G) );
    for (J=0; J<Nn; J++) TRY( MatDestroy(&mats_G[J]) );
    TRY( PetscFree(mats_G) );
  }

  TRY( PetscObjectSetName((PetscObject)G,"G") );
  *G_new = G;
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTDualize"
PetscErrorCode QPTDualize(QP qp,MatInvType invType,MatRegularizationType regType)
{
  MPI_Comm         comm;
  QP               child;
  Mat              R,B,Bt,F,G,K,Kplus,Kplus_orig;
  Vec              c,d,e,f,lb,tprim,lambda;
  PetscBool        B_explicit = PETSC_FALSE, B_extension = PETSC_FALSE, B_view_spectra = PETSC_FALSE;
  PetscBool        B_nest_extension = PETSC_TRUE;
  PetscBool        mp = PETSC_FALSE;
  PetscBool        true_mp = PETSC_FALSE;

  PetscFunctionBeginI;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  TRY( PetscLogEventBegin(QPT_Dualize,qp,0,0,0) );
  TRY( QPTransformBegin(QPTDualize, QPTDualizePostSolve_Private,NULL, QP_DUPLICATE_DO_NOT_COPY,&qp,&child,&comm) );
  TRY( QPAppendOptionsPrefix(child,"dual_") );

  Kplus_orig = NULL;
  B = NULL;
  c = qp->c;
  lambda = qp->lambda;
  tprim = qp->xwork;
  K = qp->A;
  f = qp->b;
  if (!qp->B) FLLOP_SETERRQ_WORLD(PETSC_ERR_ARG_NULL,"lin. equality and/or inequality constraint matrix (BE/BI) is needed for dualization");

  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_dualize_B_explicit",&B_explicit,NULL) );
  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_dualize_B_extension",&B_extension,NULL) );
  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_dualize_B_nest_extension",&B_nest_extension,NULL) );
  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_dualize_B_view_spectra",&B_view_spectra,NULL) );

  if (B_view_spectra) {
    TRY( QPTDualizeViewBSpectra_Private(qp->B) );
  }

  TRY( MatPrintInfo(qp->B) );

  TRY( PetscLogEventBegin(QPT_Dualize_PrepareBt,qp,0,0,0) );
  if (B_extension || B_nest_extension) {
    Mat B_merged;
    MatTransposeType ttype = B_explicit ? MAT_TRANSPOSE_EXPLICIT : MAT_TRANSPOSE_CHEAPEST;
    TRY( MatCreateNestPermonVerticalMerge(comm,1,&qp->B,&B_merged) );
    TRY( PermonMatTranspose(B_merged,MAT_TRANSPOSE_EXPLICIT,&Bt) );
    TRY( MatDestroy(&B_merged) );
    if (B_extension) {
      TRY( MatConvert(Bt,MATEXTENSION,MAT_INPLACE_MATRIX,&Bt) );
    } else {
      TRY( PermonMatConvertBlocks(Bt,MATEXTENSION,MAT_INPLACE_MATRIX,&Bt) );
    }
    TRY( PermonMatTranspose(Bt,ttype,&B) );
  } else {
    TRY( PermonMatTranspose(qp->B,MAT_TRANSPOSE_CHEAPEST,&Bt) );
    if (B_explicit) {
      TRY( PermonMatTranspose(Bt,MAT_TRANSPOSE_EXPLICIT,&B) );
    } else {
      /* in this case B remains the same */
      B = qp->B;
      TRY( PetscObjectReference((PetscObject)B) );
    }
  }
  TRY( FllopPetscObjectInheritName((PetscObject)B,(PetscObject)qp->B,NULL) );
  TRY( FllopPetscObjectInheritName((PetscObject)Bt,(PetscObject)qp->B,"_T") );
  TRY( PetscLogEventEnd(  QPT_Dualize_PrepareBt,qp,0,0,0) );

  if (FllopObjectInfoEnabled) {
    TRY( PetscPrintf(comm, "B and Bt after conversion:\n") );
    TRY( MatPrintInfo(B) );
    TRY( MatPrintInfo(Bt) );
  }

  /* create stiffness matrix pseudoinverse */
  TRY( MatCreateInv(K, invType, &Kplus) );
  Kplus_orig = Kplus;
  TRY( MatPrintInfo(K) );
  TRY( FllopPetscObjectInheritName((PetscObject)Kplus,(PetscObject)K,"_plus") );
  TRY( FllopPetscObjectInheritPrefix((PetscObject)Kplus,(PetscObject)child,NULL) );

  /* get or compute stiffness matrix kernel (R) */
  R = NULL;
  TRY( QPGetOperatorNullSpace(qp,&R) );
  if (R) {
    TRY( MatInvSetNullSpace(Kplus,R) );
  } else {
    TRY( PetscInfo(qp,"null space matrix not set => using -qpt_dualize_Kplus_left and -regularize 0\n") );
    mp = PETSC_TRUE;
    true_mp = PETSC_TRUE;
    regType = MAT_REG_NONE;
    TRY( PetscInfo(qp,"null space matrix not set => trying to compute one\n") );
    TRY( MatInvComputeNullSpace(Kplus) );
    TRY( MatInvGetNullSpace(Kplus,&R) );
    TRY( MatOrthColumns(R, MAT_ORTH_GS, MAT_ORTH_FORM_EXPLICIT, &R, NULL) );
    TRY( QPSetOperatorNullSpace(qp,R) );
  }
  TRY( MatInvSetRegularizationType(Kplus,regType) );
  TRY( MatSetFromOptions(Kplus) );

  TRY( PetscLogEventBegin(QPT_Dualize_FactorK,qp,Kplus,0,0) );
  TRY( MatInvSetUp(Kplus) );
  TRY( PetscLogEventEnd  (QPT_Dualize_FactorK,qp,Kplus,0,0) );

  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_dualize_Kplus_mp",&true_mp,NULL) );
  if (!true_mp) {
    TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_dualize_Kplus_left",&mp,NULL) );
  }

  /* convert to Moore-Penrose pseudoinverse using projector to image of K (kernel of R') */
  if (true_mp || mp) {
    QPPF pf_R;
    Mat P_R;
    Mat mats[3];
    Mat Kplus_new;
    Mat Rt;
    PetscInt size=2;

    if (true_mp) {
      TRY( PetscInfo(qp,"creating Moore-Penrose inverse\n") );
    } else {
      TRY( PetscInfo(qp,"creating left generalized inverse\n") );
    }
    TRY( PermonMatTranspose(R,MAT_TRANSPOSE_CHEAPEST,&Rt) );
    TRY( PetscObjectSetName((PetscObject)Rt,"Rt") );
    TRY( QPPFCreate(comm,&pf_R) );
    TRY( PetscObjectSetOptionsPrefix((PetscObject)pf_R,"Kplus_") );
    TRY( QPPFSetG(pf_R,Rt) );
    TRY( QPPFCreateP(pf_R,&P_R) );
    TRY( QPPFSetUp(pf_R) );

    if (true_mp) {
      size=3;
      mats[2]=P_R;
    }
    mats[1]=Kplus; mats[0]=P_R;
    TRY( MatCreateProd(comm,size,mats,&Kplus_new) );
    TRY( PetscObjectCompose((PetscObject)Kplus_new,"Kplus",(PetscObject)Kplus) );

    TRY( MatDestroy(&Kplus) );
    TRY( MatDestroy(&Rt) );
    TRY( MatDestroy(&P_R) );
    TRY( QPPFDestroy(&pf_R) );
    Kplus = Kplus_new;

    if (FllopDebugEnabled) {
      /* is Kplus MP? */
      Mat mats2[3];
      Mat prod;
      PetscBool flg;
      mats2[2]=K; mats2[1]=Kplus; mats2[0]=K;
      TRY( MatCreateProd(comm,3,mats2,&prod) );
      TRY( MatMultEqual(prod,K,3,&flg) );
      FLLOP_ASSERT(flg,"Kplus is left generalized inverse");
      TRY( MatDestroy(&prod) );
      if (true_mp) {
        mats2[2]=Kplus; mats2[1]=K; mats2[0]=Kplus;
        TRY( MatCreateProd(comm,3,mats2,&prod) );
        TRY( MatMultEqual(prod,Kplus,3,&flg) );
        FLLOP_ASSERT(flg,"Kplus is Moore-Penrose pseudoinverse");
        TRY( MatDestroy(&prod) );
      }
    }
  }

  TRY( PetscObjectSetName((PetscObject)Kplus,"Kplus") );

  /* G = R'*B' */
  TRY( PetscLogEventBegin(QPT_Dualize_AssembleG,qp,0,0,0) );
  TRY( MatTransposeMatMult_R_Bt(R,Bt,&G) );
  TRY( PetscLogEventEnd(  QPT_Dualize_AssembleG,qp,0,0,0) );

  /* F = B*Kplus*Bt (implicitly) */
  {
    Mat F_arr[3];
    
    TRY( MatCreateTimer(B,&F_arr[2]) );
    TRY( MatCreateTimer(Kplus,&F_arr[1]) );
    TRY( MatCreateTimer(Bt,&F_arr[0]) );
  
    TRY( MatCreateProd(comm, 3, F_arr, &F) );
    TRY( PetscObjectSetName((PetscObject) F, "F") );
    TRY( MatCreateTimer(F,&F) );
    
    TRY( PetscObjectCompose((PetscObject)F,"B",(PetscObject)B) );
    TRY( PetscObjectCompose((PetscObject)F,"K",(PetscObject)K) );
    TRY( PetscObjectCompose((PetscObject)F,"Kplus",(PetscObject)Kplus) );
    TRY( PetscObjectCompose((PetscObject)F,"Kplus_orig",(PetscObject)Kplus_orig) );
    TRY( PetscObjectCompose((PetscObject)F,"Bt",(PetscObject)Bt) );
    
    TRY( MatDestroy(&B) );       B     = F_arr[2];
    TRY( MatDestroy(&Kplus) );   Kplus = F_arr[1];
    TRY( MatDestroy(&Bt) );      Bt    = F_arr[0];
  }

  /* d = B*Kplus*f - c */
  TRY( VecDuplicate(lambda,&d) );
  TRY( MatMult(Kplus, f, tprim) );
  TRY( MatMult(B, tprim, d) );
  if(c) TRY( VecAXPY(d,-1.0,c) );

  /* e = R'*f */
  TRY( MatCreateVecs(R,&e,NULL) );
  TRY( MatMultTranspose(R,f,e) );

  /* lb(E) = -inf; lb(I) = 0 */
  if (qp->BI) {
    TRY( VecDuplicate(lambda,&lb) );
    if (qp->BE) {
      IS rows[2];
      Vec lbE,lbI;

      TRY( MatNestGetISs(qp->B,rows,NULL) );

      /* lb(E) = -inf */
      TRY( VecGetSubVector(lb,rows[0],&lbE) );
      TRY( VecSet(lbE,PETSC_NINFINITY) );
      TRY( VecRestoreSubVector(lb,rows[0],&lbE) );

      /* lb(I) = 0 */
      TRY( VecGetSubVector(lb,rows[1],&lbI) );
      TRY( VecSet(lbI,0.0) );
      TRY( VecRestoreSubVector(lb,rows[1],&lbI) );

      TRY( PetscObjectCompose((PetscObject)lb,"is_finite",(PetscObject)rows[1]) );
    } else {
      /* lb = 0 */
      TRY( VecSet(lb,0.0) );
    }
  } else {
    lb = NULL;
  }

  //TODO what if initial x is nonzero
  /* lambda = o */
  TRY( VecZeroEntries(lambda) );

  /* set data of the new aux. QP */
  TRY( QPSetOperator(child, F) );
  TRY( QPSetRhs(child, d) );
  TRY( QPSetEq(child, G, e) );
  TRY( QPSetIneq(child, NULL, NULL) );
  TRY( QPSetBox(child, NULL, lb, NULL) );
  TRY( QPSetInitialVector(child,lambda) );
  
  /* create special preconditioner for dual formulation */
  {
    TRY( PCDestroy(&child->pc) );
    TRY( QPGetPC(child,&child->pc) );
    TRY( PCSetType(child->pc,PCDUAL) );
    TRY( FllopPetscObjectInheritPrefix((PetscObject)child->pc,(PetscObject)child,NULL) );
  }

  TRY( PetscLogEventEnd(QPT_Dualize,qp,0,0,0) );

  TRY( QPTDualizeView_Private(qp,child) );

  TRY( MatDestroy(&B) );
  TRY( MatDestroy(&Bt) );
  TRY( MatDestroy(&F) );
  TRY( MatDestroy(&G) );
  TRY( MatDestroy(&Kplus) );
  TRY( VecDestroy(&d) );
  TRY( VecDestroy(&e) );
  TRY( VecDestroy(&lb) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTFetiPrepare"
PetscErrorCode QPTFetiPrepare(QP qp,PetscBool regularize)
{
  PetscFunctionBeginI;
  TRY( PetscLogEventBegin(QPT_FetiPrepare,qp,0,0,0) );
  TRY( QPTDualize(qp, MAT_INV_BLOCKDIAG, regularize ? MAT_REG_EXPLICIT : MAT_REG_NONE) );
  TRY( QPTHomogenizeEq(qp) );
  TRY( QPTEnforceEqByProjector(qp) );
  TRY( PetscLogEventEnd  (QPT_FetiPrepare,qp,0,0,0) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTFetiPrepareReuseCP"
PetscErrorCode QPTFetiPrepareReuseCP(QP qp,PetscBool regularize)
{
  QP dualQP;

  PetscFunctionBeginI;
  TRY( QPTDualize(qp, MAT_INV_BLOCKDIAG, regularize ? MAT_REG_EXPLICIT : MAT_REG_NONE) );

  TRY( QPChainGetLast(qp, &dualQP) );
  if (QPReusedCP) {
    Mat G;
    Vec e;

    /* reuse the coarse problem from the 0th iteration */
    TRY( QPSetQPPF(dualQP, QPReusedCP) );
    TRY( QPPFGetG(QPReusedCP, &G) );
    TRY( QPGetEq(dualQP, NULL, &e) );
    TRY( QPSetEq(dualQP, G, e) );
  } else {
    /* store the coarse problem from the 0th iteration */
    TRY( QPGetQPPF(dualQP, &QPReusedCP) );
    TRY( QPPFSetUp(QPReusedCP) );
    TRY( PetscObjectReference((PetscObject) QPReusedCP) );
  }

  TRY( QPTHomogenizeEq(qp) );
  TRY( QPTEnforceEqByProjector(qp) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTFetiPrepareReuseCPReset"
PetscErrorCode QPTFetiPrepareReuseCPReset()
{
  PetscFunctionBegin;
  TRY( QPPFDestroy(&QPReusedCP) );
  QPReusedCP = NULL;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolve_QPTRemoveGluingOfDirichletDofs"
static PetscErrorCode QPTPostSolve_QPTRemoveGluingOfDirichletDofs(QP child,QP parent)
{
  IS is = (IS) child->postSolveCtx;
  IS *iss_parent, *iss_child;
  Vec *lambda_parent, *lambda_child;
  PetscInt Mn,II;

  PetscFunctionBegin;
  Mn = child->BE_nest_count;
  FLLOP_ASSERT(Mn>=1,"child->BE_nest_count >= 2");
  TRY( PetscMalloc4(Mn,&iss_parent,Mn,&iss_child,Mn,&lambda_parent,Mn,&lambda_child) );

  TRY( MatNestGetISs(parent->BE,iss_parent,NULL) );
  TRY( MatNestGetISs(child->BE,iss_child,NULL) );

  /* 0,1 corresponds to gluing and Dirichlet part, respectively */
  for (II=0; II<Mn; II++) {
    TRY( VecGetSubVector(parent->lambda_E,iss_parent[II],&lambda_parent[II]) );
    TRY( VecGetSubVector(child->lambda_E, iss_child[II], &lambda_child[II]) );
  }

#if 0
  {
    PetscInt start, end;
    IS is_removed;

    TRY( VecGetOwnershipRange(lambda_parent[0],&start,&end) );
    TRY( ISComplement(is,start,end,&is_removed) );

    TRY( PetscPrintf(PETSC_COMM_WORLD, "#### "__FUNCT__": is_removed:\n") );
    TRY( ISView(is_removed,PETSC_VIEWER_STDOUT_WORLD) );

    TRY( PetscPrintf(PETSC_COMM_WORLD, "\n\n") );
  }
#endif

  /* copy values from the multiplier related to the restricted Bg to the multiplier
     corresponding to the original Bg, pad with zeros */
  TRY( VecZeroEntries(lambda_parent[0]) );
  //TODO PETSc 3.5+ VecGetSubVector,VecRestoreSubVector
  {
    VecScatter sc;
    TRY( VecScatterCreate(lambda_child[0],NULL,lambda_parent[0],is,&sc) );
    TRY( VecScatterBegin(sc,lambda_child[0],lambda_parent[0],INSERT_VALUES,SCATTER_FORWARD) );
    TRY( VecScatterEnd(  sc,lambda_child[0],lambda_parent[0],INSERT_VALUES,SCATTER_FORWARD) );
    TRY( VecScatterDestroy(&sc) );
  }

  /* Bd is the same, just copy values */
  for (II=1; II<Mn; II++) {
    TRY( VecCopy(lambda_child[II],lambda_parent[II]) );
  }

  for (II=0; II<Mn; II++) {
    TRY( VecRestoreSubVector(parent->lambda_E,iss_parent[II],&lambda_parent[II]) );
    TRY( VecRestoreSubVector(child->lambda_E, iss_child[II], &lambda_child[II]) );
  }

  TRY( PetscFree4(iss_parent,iss_child,lambda_parent,lambda_child) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolveDestroy_QPTRemoveGluingOfDirichletDofs"
static PetscErrorCode QPTPostSolveDestroy_QPTRemoveGluingOfDirichletDofs(void *ctx)
{
  IS is = (IS) ctx;

  PetscFunctionBegin;
  TRY( ISDestroy(&is) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTRemoveGluingOfDirichletDofs"
PetscErrorCode QPTRemoveGluingOfDirichletDofs(QP qp)
{
  QP child;
  MPI_Comm comm;
  PetscBool flg;
  Mat Bg, Bgt, Bgt_new, Bg_new, Bd, Bdt;
  Mat **mats;
  PetscInt Mn,Nn,II;
  IS is;

  PetscFunctionBeginI;
  TRY( PetscObjectTypeCompareAny((PetscObject)qp->BE,&flg,MATNEST,MATNESTPERMON,"") );
  if (!flg) FLLOP_SETERRQ(PetscObjectComm((PetscObject)qp),PETSC_ERR_SUP,"only for eq. con. matrix qp->BE of type MATNEST or MATNESTPERMON");

  TRY( PetscLogEventBegin(QPT_RemoveGluingOfDirichletDofs,qp,0,0,0) );
  TRY( MatNestGetSize(qp->BE,&Mn,&Nn) );
  TRY( MatNestGetSubMats(qp->BE,&Mn,&Nn,&mats) );
  FLLOP_ASSERT(Mn>=2,"Mn==2");
  FLLOP_ASSERT(Nn==1,"Nn==1");
  FLLOP_ASSERT(!qp->cE,"!qp->cE");
  Bg = mats[0][0];
  Bd = mats[1][0];

  TRY( QPTransformBegin(QPTRemoveGluingOfDirichletDofs,
      QPTPostSolve_QPTRemoveGluingOfDirichletDofs, QPTPostSolveDestroy_QPTRemoveGluingOfDirichletDofs,
      QP_DUPLICATE_COPY_POINTERS, &qp, &child, &comm) );
  TRY( VecDestroy(&child->lambda_E) );
  TRY( MatDestroy(&child->B) );
  TRY( MatDestroy(&child->BE) );
  child->BE_nest_count = 0;

  TRY( PermonMatTranspose(Bg,MAT_TRANSPOSE_CHEAPEST,&Bgt) );
  TRY( PermonMatTranspose(Bd,MAT_TRANSPOSE_CHEAPEST,&Bdt) );

  flg = PETSC_FALSE;
  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_remove_gluing_dirichlet_old",&flg,NULL) );
  if (flg) {
    FLLOP_EXTERN PetscErrorCode MatRemoveGluingOfDirichletDofs_old(Mat,Vec,Mat,Mat*,Vec*,IS*);
    TRY( MatRemoveGluingOfDirichletDofs_old(Bgt,NULL,Bdt,&Bgt_new,NULL,&is) );
  } else {
    TRY( MatRemoveGluingOfDirichletDofs(Bgt,NULL,Bdt,&Bgt_new,NULL,&is) );
  }

  TRY( MatDestroy(&Bdt) );

  TRY( PermonMatTranspose(Bgt_new,MAT_TRANSPOSE_CHEAPEST,&Bg_new) );
  TRY( QPAddEq(child,Bg_new,NULL) );
  for (II=1; II<Mn; II++) {
    TRY( QPAddEq(child,mats[II][0],NULL) );
  }
  TRY( MatDestroy(&Bgt_new) );

  //TODO use VecGetSubVector with PETSc 3.5+
  TRY( MatCreateVecs(child->BE,NULL,&child->lambda_E) );
  TRY( VecInvalidate(child->lambda_E) );
  child->postSolveCtx = (void*) is;

  TRY( MatPrintInfo(Bg) );
  TRY( VecPrintInfo(qp->lambda_E) );
  TRY( MatPrintInfo(Bg_new) );
  TRY( VecPrintInfo(child->lambda_E) );

  TRY( MatDestroy(&Bg_new) );
  TRY( PetscLogEventEnd(QPT_RemoveGluingOfDirichletDofs,qp,0,0,0) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolve_QPTScale"
static PetscErrorCode QPTPostSolve_QPTScale(QP child,QP parent)
{
  QPTScale_Ctx *ctx = (QPTScale_Ctx*) child->postSolveCtx;

  PetscFunctionBegin;
  if (ctx->dE) {
    TRY( VecPointwiseMult(parent->lambda_E,ctx->dE,child->lambda_E) );
  }
  if (ctx->dI) {
    TRY( VecPointwiseMult(parent->lambda_I,ctx->dI,child->lambda_I) );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPTPostSolveDestroy_QPTScale"
static PetscErrorCode QPTPostSolveDestroy_QPTScale(void *ctx)
{
  QPTScale_Ctx *cctx = (QPTScale_Ctx*) ctx;  

  PetscFunctionBegin;
  TRY( VecDestroy(&cctx->dE) );
  TRY( VecDestroy(&cctx->dI) );
  TRY( VecDestroy(&cctx->dO) );
  TRY( PetscFree(cctx) );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTScale_Private"
PetscErrorCode QPTScale_Private(Mat A,Vec b,Vec d,Mat *DA,Vec *Db)
{
  PetscFunctionBegin;
  TRY( MatDuplicate(A,MAT_COPY_VALUES,DA) );
  TRY( FllopPetscObjectInheritName((PetscObject)*DA,(PetscObject)A,NULL) );
  TRY( MatDiagonalScale(*DA,d,NULL) );
  
  if (b) {
    TRY( VecDuplicate(b,Db) );
    TRY( VecPointwiseMult(*Db,d,b) );
    TRY( FllopPetscObjectInheritName((PetscObject)*Db,(PetscObject)b,NULL) );
  } else {
    *Db = NULL;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTScale"
PetscErrorCode QPTScale(QP qp)
{
  MPI_Comm comm;
  QPScaleType ScalType;
  MatOrthType R_orth_type=MAT_ORTH_GS;
  MatOrthForm R_orth_form=MAT_ORTH_FORM_EXPLICIT;
  PetscBool remove_gluing_of_dirichlet=PETSC_FALSE;
  PetscBool set;
  QP child;
  Mat A,DA;
  Vec b,d,Db;
  QPTScale_Ctx *ctx;

  PetscFunctionBeginI;
  TRY( QPChainGetLast(qp,&qp) );

  _fllop_ierr = PetscObjectOptionsBegin((PetscObject)qp);CHKERRQ(_fllop_ierr);
  TRY( PetscOptionsBool("-qp_E_remove_gluing_of_dirichlet","remove gluing of DOFs on Dirichlet boundary","QPTRemoveGluingOfDirichletDofs",remove_gluing_of_dirichlet,&remove_gluing_of_dirichlet,NULL) );
  _fllop_ierr = PetscOptionsEnd();CHKERRQ(_fllop_ierr);
  TRY( PetscInfo1(qp, "-qp_E_remove_gluing_of_dirichlet %d\n",remove_gluing_of_dirichlet) );
  if (remove_gluing_of_dirichlet) {
    TRY( QPTRemoveGluingOfDirichletDofs(qp) );
  }

  TRY( QPTransformBegin(QPTScale,
      QPTPostSolve_QPTScale, QPTPostSolveDestroy_QPTScale,
      QP_DUPLICATE_COPY_POINTERS, &qp, &child, &comm) );
  TRY( PetscNew(&ctx) );

  _fllop_ierr = PetscObjectOptionsBegin((PetscObject)qp);CHKERRQ(_fllop_ierr);
  A = qp->A;
  b = qp->b;
  d = NULL;
  ScalType = QP_SCALE_NONE;
  TRY( PetscOptionsEnum("-qp_O_scale_type", "", "QPSetSystemScaling", QPScaleTypes, (PetscEnum)ScalType, (PetscEnum*)&ScalType, &set) );
  TRY( PetscInfo1(qp, "-qp_O_scale_type %s\n",QPScaleTypes[ScalType]) );
  if (ScalType) {
    if (ScalType == QP_SCALE_ROWS_NORM_2) {
      TRY( MatGetRowNormalization(A,&d) );
    } else {
      FLLOP_SETERRQ1(comm,PETSC_ERR_SUP,"-qp_O_scale_type %s not supported",QPScaleTypes[ScalType]);
    }

    TRY( QPTScale_Private(A,b,d,&DA,&Db) );
    
    TRY( QPSetOperator(child,DA) );
    TRY( QPSetRhs(child,Db) );
    ctx->dO = d;

    TRY( MatDestroy(&DA) );
    TRY( VecDestroy(&Db) );
  }

  A = qp->BE;
  b = qp->cE;
  d = NULL;
  ScalType = QP_SCALE_NONE;
  TRY( PetscOptionsEnum("-qp_E_scale_type", "", "QPSetEqScaling", QPScaleTypes, (PetscEnum)ScalType, (PetscEnum*)&ScalType, &set) );
  TRY( PetscInfo1(qp, "-qp_E_scale_type %s\n",QPScaleTypes[ScalType]) );
  if (ScalType) {
    if (ScalType == QP_SCALE_ROWS_NORM_2) {
      TRY( MatGetRowNormalization(A,&d) );
    } else if (ScalType == QP_SCALE_DDM_MULTIPLICITY) {
      TRY( QPGetEqMultiplicityScaling(qp,&d,&ctx->dI) );
    } else {
      FLLOP_SETERRQ1(comm,PETSC_ERR_SUP,"-qp_E_scale_type %s not supported",QPScaleTypes[ScalType]);
    }

    TRY( QPTScale_Private(A,b,d,&DA,&Db) );
    
    TRY( QPSetQPPF(child,NULL) );
    TRY( QPSetEq(child,DA,Db) );
    TRY( QPSetEqMultiplier(child,NULL) );
    ctx->dE = d;

    TRY( MatDestroy(&DA) );
    TRY( VecDestroy(&Db) );
  }

  A = qp->BI;
  b = qp->cI;
  d = ctx->dI;
  ScalType = QP_SCALE_NONE;
  TRY( PetscOptionsEnum("-qp_I_scale_type", "", "QPSetIneqScaling", QPScaleTypes, (PetscEnum)ScalType, (PetscEnum*)&ScalType, &set) );
  TRY( PetscInfo1(qp, "-qp_I_scale_type %s\n",QPScaleTypes[ScalType]) );
  if (ScalType || d) {
    if (ScalType && d) FLLOP_SETERRQ1(comm,PETSC_ERR_SUP,"-qp_I_scale_type %s not supported for given eq. con. scaling",QPScaleTypes[ScalType]);

    if (ScalType == QP_SCALE_ROWS_NORM_2) {
      TRY( MatGetRowNormalization(A,&d) );
    } else if (!d) {
      FLLOP_SETERRQ1(comm,PETSC_ERR_SUP,"-qp_I_scale_type %s not supported",QPScaleTypes[ScalType]);
    }

    TRY( QPTScale_Private(A,b,d,&DA,&Db) );

    TRY( QPSetIneq(child,DA,Db) );
    TRY( QPSetIneqMultiplier(child,NULL) );
    ctx->dI = d;

    TRY( MatDestroy(&DA) );
    TRY( VecDestroy(&Db) );
  }

  if (qp->R) {
    TRY( PetscOptionsEnum("-qp_R_orth_type", "type of nullspace matrix orthonormalization", "", MatOrthTypes, (PetscEnum)R_orth_type, (PetscEnum*)&R_orth_type, NULL) );
    TRY( PetscOptionsEnum("-qp_R_orth_form", "form of nullspace matrix orthonormalization", "", MatOrthForms, (PetscEnum)R_orth_form, (PetscEnum*)&R_orth_form, NULL) );
    TRY( PetscInfo1(qp, "-qp_R_orth_type %s\n",MatOrthTypes[R_orth_type]) );
    TRY( PetscInfo1(qp, "-qp_R_orth_form %s\n",MatOrthForms[R_orth_form]) );
    if (R_orth_type) {
      Mat Rnew;
      TRY( MatOrthColumns(qp->R, R_orth_type, R_orth_form, &Rnew, NULL) );
      TRY( QPSetOperatorNullSpace(child,Rnew) );
      TRY( MatDestroy(&Rnew) );
    }
  }

  _fllop_ierr = PetscOptionsEnd();CHKERRQ(_fllop_ierr);
  child->postSolveCtx = ctx;
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTNormalizeObjective"
PetscErrorCode QPTNormalizeObjective(QP qp)
{
  PetscReal norm_A, norm_b;

  PetscFunctionBeginI;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  TRY( QPChainGetLast(qp,&qp) );
  TRY( MatGetMaxEigenvalue(qp->A, NULL, &norm_A, PETSC_DECIDE, PETSC_DECIDE) );
  TRY( VecNorm(qp->b,NORM_2,&norm_b) );
  TRY( PetscInfo2(qp,"||A||=%.12e, scale A by 1/||A||=%.12e\n",norm_A,1.0/norm_A) );
  TRY( PetscInfo2(qp,"||b||=%.12e, scale b by 1/||b||=%.12e\n",norm_b,1.0/norm_b) );
  TRY( QPTScaleObjectiveByScalar(qp, 1.0/norm_A, 1.0/norm_b) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTNormalizeHessian"
PetscErrorCode QPTNormalizeHessian(QP qp)
{
  PetscReal norm_A;

  PetscFunctionBeginI;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  TRY( QPChainGetLast(qp,&qp) );
  TRY( MatGetMaxEigenvalue(qp->A, NULL, &norm_A, PETSC_DECIDE, PETSC_DECIDE) );
  TRY( PetscInfo2(qp,"||A||=%.12e, scale A by 1/||A||=%.12e\n",norm_A,1.0/norm_A) );
  TRY( QPTScaleObjectiveByScalar(qp, 1.0/norm_A, 1.0/norm_A) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolve_QPTScaleObjectiveByScalar"
static PetscErrorCode QPTPostSolve_QPTScaleObjectiveByScalar(QP child,QP parent)
{
  QPTScaleObjectiveByScalar_Ctx *psctx = (QPTScaleObjectiveByScalar_Ctx*)child->postSolveCtx;
  PetscReal scale_A = psctx->scale_A;
  PetscReal scale_b = psctx->scale_b;
  Vec lb,ub;
  Vec llb,lub;
  Vec llbnew,lubnew;
  QPC qpc,qpcc;

  PetscFunctionBegin;
  TRY( VecCopy(child->x,parent->x) );
  TRY( VecScale(parent->x,scale_A/scale_b) );
  
  if (parent->Bt_lambda) {
    TRY( VecCopy(child->Bt_lambda,parent->Bt_lambda) );
    TRY( VecScale(parent->Bt_lambda,1.0/scale_b) );
  }
  if (parent->lambda_E) {
    TRY( VecCopy(child->lambda_E,parent->lambda_E) );
    TRY( VecScale(parent->lambda_E,1.0/scale_b) );
  }
  TRY( QPGetBox(parent,NULL,&lb,&ub) );
  TRY( QPGetQPC(parent, &qpc) );
  TRY( QPGetQPC(child, &qpcc) );
  TRY( QPCBoxGetMultipliers(qpcc,&llb,&lub) );
  TRY( QPCBoxGetMultipliers(qpc,&llbnew,&lubnew) );
  if (lb) {
    TRY( VecCopy(llb,llbnew) );
    TRY( VecScale(llbnew,1.0/scale_b) );
  }
  if (ub) {
    TRY( VecCopy(lub,lubnew) );
    TRY( VecScale(lubnew,1.0/scale_b) );
  }
  if (parent->BI) {
    TRY( VecCopy(child->lambda_I,parent->lambda_I) );
    TRY( VecScale(parent->lambda_I,1.0/scale_b) );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolveDestroy_QPTScaleObjectiveByScalar"
static PetscErrorCode QPTPostSolveDestroy_QPTScaleObjectiveByScalar(void *ctx)
{
  PetscFunctionBegin;
  PetscFree(ctx);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTScaleObjectiveByScalar"
PetscErrorCode QPTScaleObjectiveByScalar(QP qp,PetscScalar scale_A,PetscScalar scale_b)
{
  MPI_Comm comm;
  QP child;
  QPTScaleObjectiveByScalar_Ctx *ctx;
  Mat Anew;
  Vec bnew;
  Vec lb,ub;
  Vec lbnew,ubnew;
  PetscReal norm_A,norm_b;

  PetscFunctionBeginI;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  TRY( QPTransformBegin(QPTScaleObjectiveByScalar,
      QPTPostSolve_QPTScaleObjectiveByScalar, QPTPostSolveDestroy_QPTScaleObjectiveByScalar,
      QP_DUPLICATE_COPY_POINTERS, &qp, &child, &comm) );
  TRY( PetscNew(&ctx) );
  ctx->scale_A = scale_A;
  ctx->scale_b = scale_b;
  child->postSolveCtx = ctx;

  if (FllopDebugEnabled) {
    TRY( MatGetMaxEigenvalue(qp->A, NULL, &norm_A, 1e-5, 50) );
    TRY( FllopDebug1("||A||=%.12e\n",norm_A) );
    TRY( VecNorm(qp->b,NORM_2,&norm_b) );
    TRY( FllopDebug1("||b||=%.12e\n",norm_b) );
  }

  if (qp->A->ops->duplicate) {
    TRY( MatDuplicate(qp->A,MAT_COPY_VALUES,&Anew) );
  } else {
    TRY( MatCreateProd(comm,1,&qp->A,&Anew) );
  }
  TRY( MatScale(Anew,scale_A) );
  TRY( QPSetOperator(child,Anew) );
  TRY( MatDestroy(&Anew) );

  TRY( VecDuplicate(qp->b,&bnew) );
  TRY( VecCopy(qp->b,bnew) );
  TRY( VecScale(bnew,scale_b) );
  TRY( QPSetRhs(child,bnew) );
  TRY( VecDestroy(&bnew) );

  TRY( QPGetBox(qp,NULL,&lb,&ub) );
  lbnew=NULL;
  if (lb) {
    TRY( VecDuplicate(lb,&lbnew) );
    TRY( VecCopy(lb,lbnew) );
    TRY( VecScaleSkipInf(lbnew,scale_b/scale_A) );
  }
  ubnew=NULL;
  if (ub) {
    TRY( VecDuplicate(ub,&ubnew) );
    TRY( VecCopy(ub,ubnew) );
    TRY( VecScaleSkipInf(ubnew,scale_b/scale_A) );
  }
  TRY( QPSetBox(child,NULL,lbnew,ubnew) );
  TRY( VecDestroy(&lbnew) );
  TRY( VecDestroy(&ubnew) );

  TRY( QPSetInitialVector(child,NULL) );
  TRY( QPSetEqMultiplier(child,NULL) );
  TRY( QPSetIneqMultiplier(child,NULL) );

  if (FllopDebugEnabled) {
    TRY( MatGetMaxEigenvalue(child->A, NULL, &norm_A, 1e-5, 50) );
    TRY( FllopDebug1("||A_new||=%.12e\n",norm_A) );
    TRY( VecNorm(child->b,NORM_2,&norm_b) );
    TRY( FllopDebug1("||b_new||=%.12e\n",norm_b) );
  }
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolve_QPTFreezeIneq"
static PetscErrorCode QPTPostSolve_QPTFreezeIneq(QP child,QP parent)
{
  PetscInt Mn;
  IS *iss=NULL,is=NULL;
  Vec lambda_EI;

  PetscFunctionBegin;
  Mn = child->BE_nest_count;
  FLLOP_ASSERT(Mn>=1,"child->BE_nest_count >= 1");

  if (Mn>1) {
    TRY( PetscMalloc1(Mn,&iss) );
    TRY( MatNestGetISs(child->BE,iss,NULL) );
    is = iss[Mn-1];

    /* copy the corresponding part of child's lambda_E to parent's lambda_I */
    TRY( VecGetSubVector(child->lambda_E,is,&lambda_EI) );
    TRY( VecCopy(lambda_EI,parent->lambda_I) );
    TRY( VecRestoreSubVector(child->lambda_E,is,&lambda_EI) );

    /* copy the rest of child's lambda_E to parent's lambda_E */
    TRY( ISConcatenate(PetscObjectComm((PetscObject)child->BE),Mn-1,iss,&is) );
    TRY( VecGetSubVector(child->lambda_E,is,&lambda_EI) );
    TRY( VecCopy(lambda_EI,parent->lambda_E) );
    TRY( VecRestoreSubVector(child->lambda_E,is,&lambda_EI) );

    TRY( ISDestroy(&is) );
    TRY( PetscFree(iss) );
  } else { /* n==1 => there had been only ineq. con. before */
    TRY( VecCopy(child->lambda_E,parent->lambda_I) );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTFreezeIneq"
PetscErrorCode QPTFreezeIneq(QP qp)
{
  MPI_Comm          comm;
  QP                child;
  Mat               BI;
  Vec               cI;

  PetscFunctionBeginI;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  TRY( QPTransformBegin(QPTFreezeIneq, QPTPostSolve_QPTFreezeIneq,NULL, QP_DUPLICATE_COPY_POINTERS,&qp,&child,&comm) );
  TRY( QPGetIneq(qp,&BI,&cI) );
  TRY( QPAddEq(child,BI,cI) );
  TRY( QPSetEqMultiplier(child,NULL) );
  TRY( QPSetIneq(child,NULL,NULL) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTSplitBE"
PetscErrorCode QPTSplitBE(QP qp)
{
  QP child;
  MPI_Comm comm;
  Mat Be=NULL, Bd=NULL, Bg=NULL;
  Mat Bet=NULL, Bdt=NULL, Bgt=NULL;
  PetscInt i, ilo, ihi, j, k, ncols, ng=0, nd=0;
  PetscInt *idxg, *idxd;
  const PetscInt *cols;
  const PetscScalar *vals;
  IS isrowg, isrowd;
  PetscBool flg;
  
  PetscFunctionBeginI;
  TRY( PetscObjectGetComm((PetscObject)qp,&comm) );
  FLLOP_ASSERT(!qp->cE,"!qp->cE");
  TRY( MatIsImplicitTranspose(qp->BE, &flg) );
  FLLOP_ASSERT(flg,"BE is implicit transpose");
  
  TRY( PetscLogEventBegin(QPT_SplitBE,qp,0,0,0) );
  TRY( QPTransformBegin(QPTSplitBE, NULL, NULL, QP_DUPLICATE_COPY_POINTERS, &qp, &child, &comm) );

  TRY( PermonMatTranspose(child->BE, MAT_TRANSPOSE_CHEAPEST, &Bet) );
  TRY( PermonMatTranspose(Bet, MAT_TRANSPOSE_EXPLICIT, &Be) );
  TRY( MatDestroy(&Bet) );

  TRY( MatGetOwnershipRange(Be, &ilo, &ihi) );

  TRY( PetscMalloc((ihi-ilo)*sizeof(PetscInt), &idxg) );
  TRY( PetscMalloc((ihi-ilo)*sizeof(PetscInt), &idxd) );

  for (i=ilo; i<ihi; i++){
    TRY( MatGetRow(Be, i, &ncols, &cols, &vals) );
    k = 0;
    for (j=0; j<ncols; j++){
      if (vals[j]) k++;
    }
    if (k==1){
      idxd[nd] = i;
      nd += 1;
    }else if (k>1){
      idxg[ng] = i;
      ng += 1;
    }else{
      FLLOP_SETERRQ(comm, PETSC_ERR_COR, "B columns can't have every element zero");
    }
    TRY( MatRestoreRow(Be, i, &ncols, &cols, &vals) );
  }
  TRY( ISCreateGeneral(comm, ng, idxg, PETSC_OWN_POINTER, &isrowg) );
  TRY( ISCreateGeneral(comm, nd, idxd, PETSC_OWN_POINTER, &isrowd) );
  
  TRY( MatCreateSubMatrix(Be, isrowg, NULL, MAT_INITIAL_MATRIX, &Bg) );
  TRY( MatCreateSubMatrix(Be, isrowd, NULL, MAT_INITIAL_MATRIX, &Bd) );
  TRY( MatDestroy(&Be) );

  TRY( PermonMatTranspose(Bg, MAT_TRANSPOSE_EXPLICIT, &Bgt) );
  TRY( MatDestroy(&Bg) );
  TRY( PermonMatTranspose(Bgt, MAT_TRANSPOSE_CHEAPEST, &Bg) );
  
  TRY( PermonMatTranspose(Bd, MAT_TRANSPOSE_EXPLICIT, &Bdt) );
  TRY( MatDestroy(&Bd) );
  TRY( PermonMatTranspose(Bdt, MAT_TRANSPOSE_CHEAPEST, &Bd) );

  TRY( MatDestroy(&child->BE) );
  TRY( QPAddEq(child, Bg, NULL) );
  TRY( QPAddEq(child, Bd, NULL) );
  
  TRY( PetscLogEventEnd(QPT_SplitBE,qp,0,0,0) );

  TRY( ISDestroy(&isrowg) );
  TRY( ISDestroy(&isrowd) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolve_QPTMatISToBlockDiag"
static PetscErrorCode QPTPostSolve_QPTMatISToBlockDiag(QP child,QP parent)
{
  Mat AsubCopy;
  Vec dir,b_adjust,b_adjustU,resid;
  Mat_IS *matis  = (Mat_IS*)parent->A->data;
  QPTMatISToBlockDiag_Ctx *ctx = (QPTMatISToBlockDiag_Ctx*)child->postSolveCtx;
  PetscReal norm=0.0,normb=0.0;
  PetscBool computeNorm = PETSC_FALSE;
  PetscInt lock;

  PetscFunctionBegin;
  TRY( PetscOptionsGetBool(NULL,NULL,"-qpt_matis_to_diag_norm",&computeNorm,NULL) );
  if (ctx->isDir) {
    if (computeNorm) {
      /* TODO: change implementation for submat copies for PETSc>=3.8 */
      MatDuplicate(child->A,MAT_COPY_VALUES,&AsubCopy);
      //TRY( MatGetLocalToGlobalMapping(child->A,&l2g,NULL) );
      //TRY( ISGlobalToLocalMappingApplyIS(l2g,IS_GTOLM_DROP,ctx->isDir,&isDirLoc) );
      //TRY( MatGetLocalSubMatrix(parent->A,isDirLoc,isDirLoc,&Asub) );
      //TRY( MatDuplicate(Asub,MAT_COPY_VALUES,&AsubCopy) );
      TRY( VecDuplicate(child->x,&dir) );
      TRY( VecDuplicate(child->b,&b_adjustU) );
      TRY( VecCopy(child->b,b_adjustU) );
      TRY( VecGetLocalVector(dir,matis->y) );
      TRY( VecScatterBegin(matis->cctx,parent->x,matis->y,INSERT_VALUES,SCATTER_FORWARD) ); /* set local vec */
      TRY( VecScatterEnd(matis->cctx,parent->x,matis->y,INSERT_VALUES,SCATTER_FORWARD) );
      TRY( VecRestoreLocalVector(dir,matis->y) );
      TRY( MatZeroRowsColumnsIS(child->A,ctx->isDir,1.0,dir,b_adjustU) );
    }
  } else {
    /* propagate changed RHS */
    /* TODO: flag for pure Neumann? */
    TRY( VecGetLocalVector(child->b,matis->y) );
    TRY( VecLockGet(parent->b,&lock) );
    if (lock) TRY( VecLockReadPop(parent->b) ); /* TODO: safe? */
    TRY( VecSet(parent->b,0.0) );
    TRY( VecScatterBegin(matis->rctx,matis->y,parent->b,ADD_VALUES,SCATTER_REVERSE) );
    TRY( VecScatterEnd(matis->rctx,matis->y,parent->b,ADD_VALUES,SCATTER_REVERSE) );
    if (lock) TRY( VecLockReadPush(parent->b) );
    TRY( VecRestoreLocalVector(child->b,matis->y) );
  }

  /* assemble solution */
  TRY( VecGetLocalVector(child->x,matis->x) );
  TRY( VecScatterBegin(matis->rctx,matis->x,parent->x,INSERT_VALUES,SCATTER_REVERSE) );
  TRY( VecScatterEnd(matis->rctx,matis->x,parent->x,INSERT_VALUES,SCATTER_REVERSE) );
  TRY( VecRestoreLocalVector(child->x,matis->x) );
  
  if (computeNorm) {
    /* compute norm */
    TRY( VecDuplicate(parent->b,&resid) );
    if (ctx->isDir) {
      TRY( VecDuplicate(parent->b,&b_adjust) );
      TRY( VecSet(b_adjust,.0) );
      TRY( VecGetLocalVector(b_adjustU,matis->y) );
      TRY( VecScatterBegin(matis->rctx,matis->y,b_adjust,ADD_VALUES,SCATTER_REVERSE) );
      TRY( VecScatterEnd(matis->rctx,matis->y,b_adjust,ADD_VALUES,SCATTER_REVERSE) );
      TRY( VecRestoreLocalVector(b_adjustU,matis->y) );
      TRY( MatMult(parent->A,parent->x,resid) );
      TRY( VecAXPY(resid,-1.0,b_adjust) ); /* Ax-b */ 
      TRY( VecNorm(b_adjust,NORM_2,&normb) );
      TRY( MatCopy(AsubCopy,child->A,SAME_NONZERO_PATTERN) );
      //TRY( MatCopy(AsubCopy,Asub,SAME_NONZERO_PATTERN) );
      TRY( VecDestroy(&b_adjustU) );
      TRY( VecDestroy(&b_adjust) );
    } else {
      TRY( MatMult(parent->A,parent->x,resid) );
      TRY( VecAXPY(resid,-1.0,parent->b) ); /* Ax-b */ 
      TRY( VecNorm(parent->b,NORM_2,&normb) );
    } 
    TRY( VecNorm(resid,NORM_2,&norm) );
    TRY( PetscPrintf(PETSC_COMM_WORLD,"Dirichlet in Hess: %d, r = ||Ax-b|| = %e, r/||b|| = %e\n",!ctx->isDir,norm,norm/normb) );
    TRY( VecDestroy(&resid) );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTPostSolveDestroy_QPTMatISToBlockDiag"
static PetscErrorCode QPTPostSolveDestroy_QPTMatISToBlockDiag(void *ctx)
{
  QPTMatISToBlockDiag_Ctx *cctx = (QPTMatISToBlockDiag_Ctx*)ctx;

  PetscFunctionBegin;
  TRY( ISDestroy(&cctx->isDir) );
  TRY( PetscFree(cctx) );
  PetscFunctionReturn(0);
}

/*@
   QPTMatISToBlockDiag - Transforms system matrix from MATIS format to BlockDiag.

   Collective on QP

   Input Parameter:
.  qp   - the QP

   Level: developer
@*/
#undef __FUNCT__
#define __FUNCT__ "QPTMatISToBlockDiag"
PetscErrorCode QPTMatISToBlockDiag(QP qp)
{
  QP child;
  Mat A;
  IS l2g,i2g;
  const PetscInt *idx_l2g;
  QPTMatISToBlockDiag_Ctx *ctx;
  ISLocalToGlobalMapping mapping;
  PetscInt n;         /* number of nodes (interior+interface) in this subdomain */
  PetscInt n_B;       /* number of interface nodes in this subdomain */
  PetscInt n_I;
  PetscInt n_neigh;   /* number of neighbours this subdomain has (by now, INCLUDING OR NOT the subdomain itself). */
                      /* Once this is definitively decided, the code can be simplifies and some if's eliminated.  */
  PetscInt *neigh;    /* list of neighbouring subdomains                                                          */
  PetscInt *n_shared; /* n_shared[j] is the number of nodes shared with subdomain neigh[j]                        */
  PetscInt **shared;  /* shared[j][i] is the local index of the i-th node shared with subdomain neigh[j]          */
  PetscInt *idx_I_local,*idx_B_local,*idx_I_global,*idx_B_global;
  PetscInt i,j;
  PetscBT  bt;
  IS is_B_local,is_I_local,is_B_global, is_I_global; /* local (seq) index sets for interface (B) and interior (I) nodes */
  VecScatter N_to_B;      /* scattering context from all local nodes to local interface nodes */
  VecScatter global_to_B; /* scattering context from global to local interface nodes */
  Vec b,D,vec1_B;
  MPI_Comm comm;
  
  PetscFunctionBeginI;
  PetscValidHeaderSpecific(qp,QP_CLASSID,1);
  TRY( PetscObjectGetComm((PetscObject)qp,&comm) );
  TRY( QPTransformBegin(QPTMatISToBlockDiag,QPTPostSolve_QPTMatISToBlockDiag,QPTPostSolveDestroy_QPTMatISToBlockDiag,QP_DUPLICATE_DO_NOT_COPY,&qp,&child,&comm) );
  TRY( PetscNew(&ctx) );
  child->postSolveCtx = ctx;
  ctx->isDir = NULL; /* set by QPFetiSetDirichlet() */
  TRY( PCDestroy(&child->pc) );
  TRY( QPGetPC(child,&child->pc) );

  /* create block diag */
  Mat_IS *matis  = (Mat_IS*)qp->A->data;
  TRY( MatCreateBlockDiag(comm,matis->A,&A) );
  TRY( QPSetOperator(child,A) );
  TRY( QPSetEq(child,qp->BE,NULL) );
  TRY( QPSetOperatorNullSpace(child,qp->R) );
  TRY( MatDestroy(&qp->BE) );

  /* get mappings for RHS decomposition
  *  create interface mappings
  *  adapted from PCISSetUp */
  /* get info on mapping */
  mapping = qp->A->rmap->mapping;
  TRY( ISLocalToGlobalMappingGetSize(mapping,&n) );
  TRY( ISLocalToGlobalMappingGetInfo(mapping,&n_neigh,&neigh,&n_shared,&shared) );

  /* Identifying interior and interface nodes, in local numbering */
  TRY( PetscBTCreate(n,&bt) );
  for (i=0;i<n_neigh;i++) {
    for (j=0;j<n_shared[i];j++) {
        TRY( PetscBTSet(bt,shared[i][j]) );
    }
  }
  /* Creating local and global index sets for interior and inteface nodes. */
  TRY( PetscMalloc1(n,&idx_I_local) );
  TRY( PetscMalloc1(n,&idx_B_local) );
  for (i=0, n_B=0, n_I=0; i<n; i++) {
    if (!PetscBTLookup(bt,i)) {
      idx_I_local[n_I] = i;
      n_I++;
    } else {
      idx_B_local[n_B] = i;
      n_B++;
    }
  }
  /* Getting the global numbering */
  idx_B_global = idx_I_local + n_I; /* Just avoiding allocating extra memory, since we have vacant space */
  idx_I_global = idx_B_local + n_B;
  TRY( ISLocalToGlobalMappingApply(mapping,n_B,idx_B_local,idx_B_global) );
  TRY( ISLocalToGlobalMappingApply(mapping,n_I,idx_I_local,idx_I_global) );

  /* Creating the index sets */
  TRY( ISCreateGeneral(PETSC_COMM_SELF,n_B,idx_B_local,PETSC_COPY_VALUES, &is_B_local) );
  TRY( ISCreateGeneral(PETSC_COMM_SELF,n_B,idx_B_global,PETSC_COPY_VALUES,&is_B_global) );
  /* TODO remove interior idx sets */
  TRY( ISCreateGeneral(PETSC_COMM_SELF,n_I,idx_I_local,PETSC_COPY_VALUES, &is_I_local) );
  TRY( ISCreateGeneral(PETSC_COMM_SELF,n_I,idx_I_global,PETSC_COPY_VALUES,&is_I_global) );

    /* Creating work vectors and arrays */
  TRY( VecCreateSeq(PETSC_COMM_SELF,n_B,&vec1_B) );
  TRY( VecDuplicate(vec1_B, &D) );

  /* Creating the scatter contexts */
  TRY( VecScatterCreate(matis->x,is_B_local,vec1_B,NULL,&N_to_B) );
  TRY( VecScatterCreate(qp->x,is_B_global,vec1_B,NULL,&global_to_B) );

  /* Creating scaling "matrix" D */
  TRY( VecSet(D,1.0) );
  TRY( VecScatterBegin(N_to_B,matis->counter,vec1_B,INSERT_VALUES,SCATTER_FORWARD) );
  TRY( VecScatterEnd(N_to_B,matis->counter,vec1_B,INSERT_VALUES,SCATTER_FORWARD) );
  TRY( VecPointwiseDivide(D,D,vec1_B) );

  /* decompose assembled vecs */
  TRY( MatCreateVecs(A,&child->x,&child->b) );
  /* assemble b */
  TRY( VecGetLocalVector(child->b,matis->y) );
  TRY( VecDuplicate(qp->b,&b) );
  TRY( VecCopy(qp->b,b) );
  TRY( VecScatterBegin(global_to_B,qp->b,vec1_B,INSERT_VALUES,SCATTER_FORWARD) ); /* get interface DOFs */
  TRY( VecScatterEnd(global_to_B,qp->b,vec1_B,INSERT_VALUES,SCATTER_FORWARD) );
  TRY( VecPointwiseMult(vec1_B,D,vec1_B) ); /* DOF/(number of subdomains it belongs to) */
  TRY( VecScatterBegin(global_to_B,vec1_B,b,INSERT_VALUES,SCATTER_REVERSE) ); /* replace values in RHS */
  TRY( VecScatterEnd(global_to_B,vec1_B,b,INSERT_VALUES,SCATTER_REVERSE) );
  TRY( VecScatterBegin(matis->cctx,b,matis->y,INSERT_VALUES,SCATTER_FORWARD) ); /* set local vec */
  TRY( VecScatterEnd(matis->cctx,b,matis->y,INSERT_VALUES,SCATTER_FORWARD) );
  /* assemble x */
  TRY( VecGetLocalVector(child->x,matis->x) );
  TRY( VecScatterBegin(matis->cctx,qp->x,matis->x,INSERT_VALUES,SCATTER_FORWARD) ); /* set local vec */
  TRY( VecScatterEnd(matis->cctx,qp->x,matis->x,INSERT_VALUES,SCATTER_FORWARD) );

  /* inherit l2g and i2g */
  TRY( ISLocalToGlobalMappingGetIndices(mapping,&idx_l2g) );
  TRY( ISCreateGeneral(PetscObjectComm((PetscObject)qp),n,idx_l2g,PETSC_COPY_VALUES,&l2g) );
  TRY( ISLocalToGlobalMappingRestoreIndices(mapping,&idx_l2g) );
  TRY( QPFetiSetLocalToGlobalMapping(child,l2g) );
  TRY( ISOnComm(is_B_global,PETSC_COMM_WORLD,PETSC_COPY_VALUES,&i2g) );
  TRY( ISSort(i2g) );
  TRY( QPFetiSetInterfaceToGlobalMapping(child,i2g) );

  TRY( ISDestroy(&i2g) );
  TRY( ISDestroy(&l2g) );
  TRY( VecRestoreLocalVector(child->x,matis->x) );
  TRY( VecRestoreLocalVector(child->b,matis->y) );
  TRY( ISLocalToGlobalMappingRestoreInfo(mapping,&n_neigh,&neigh,&n_shared,&shared) );
  TRY( ISDestroy(&is_B_local) );
  TRY( ISDestroy(&is_B_global) );
  TRY( ISDestroy(&is_I_local) );
  TRY( ISDestroy(&is_I_global) );
  TRY( VecScatterDestroy(&N_to_B) );
  TRY( VecScatterDestroy(&global_to_B) );
  TRY( PetscFree(idx_B_local) );
  TRY( PetscFree(idx_I_local) );
  TRY( PetscBTDestroy(&bt) );
  TRY( VecDestroy(&D) );
  TRY( VecDestroy(&vec1_B) );
  TRY( VecDestroy(&b) );
  TRY( MatDestroy(&A) );
  
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTAllInOne"
PetscErrorCode QPTAllInOne(QP qp,MatInvType invType,PetscBool dual,PetscBool project,PetscReal penalty,PetscBool penalty_direct,PetscBool regularize)
{
  MatRegularizationType regularize_e;
  PetscBool freeze=PETSC_FALSE, normalize=PETSC_FALSE, normalize_hessian=PETSC_FALSE;
  QP last;

  PetscFunctionBeginI;
  regularize_e = regularize ? MAT_REG_EXPLICIT : MAT_REG_NONE;

  TRY( PetscLogEventBegin(QPT_AllInOne,qp,0,0,0) );
  _fllop_ierr = PetscObjectOptionsBegin((PetscObject)qp);CHKERRQ(_fllop_ierr);
  TRY( PetscOptionsBool("-qp_I_freeze","perform QPTFreezeIneq","QPTFreezeIneq",freeze,&freeze,NULL) );
  TRY( PetscOptionsBoolGroupBegin("-qp_O_normalize","perform QPTNormalizeObjective","QPTNormalizeObjective",&normalize) );
  TRY( PetscOptionsBoolGroupEnd("-qp_O_normalize_hessian","perform QPTNormalizeHessian","QPTNormalizeHessian",&normalize_hessian) );
  _fllop_ierr = PetscOptionsEnd();CHKERRQ(_fllop_ierr);

  //TODO do this until QPTFromOptions supports chain updates
  TRY( QPRemoveChild(qp) );

  if (normalize) {
    TRY( QPTNormalizeObjective(qp) );
  } else if (normalize_hessian) {
    TRY( QPTNormalizeHessian(qp) );
  }

  TRY( QPTScale(qp) );
  TRY( QPTOrthonormalizeEqFromOptions(qp) );
  if (freeze) TRY( QPTFreezeIneq(qp) );
  if (dual) {
    TRY( QPTDualize(qp,invType,regularize_e) );
    TRY( QPTScale(qp) );
    TRY( QPTOrthonormalizeEqFromOptions(qp) );
  }
  if (project) {
    TRY( QPTEnforceEqByProjector(qp) );
  }
  
  if (dual || project) {
    normalize = PETSC_FALSE;
    normalize_hessian = PETSC_FALSE;
    TRY( QPChainGetLast(qp,&last) );
    _fllop_ierr = PetscObjectOptionsBegin((PetscObject)last);CHKERRQ(_fllop_ierr);
    TRY( PetscOptionsBoolGroupBegin("-qp_O_normalize","perform QPTNormalizeObjective","QPTNormalizeObjective",&normalize) );
    TRY( PetscOptionsBoolGroupEnd("-qp_O_normalize_hessian","perform QPTNormalizeHessian","QPTNormalizeHessian",&normalize_hessian) );
    _fllop_ierr = PetscOptionsEnd();CHKERRQ(_fllop_ierr);
    if (normalize) {
      TRY( QPTNormalizeObjective(qp) );
    } else if (normalize_hessian) {
      TRY( QPTNormalizeHessian(qp) );
    }
  }

  TRY( QPTEnforceEqByPenalty(qp,penalty,penalty_direct) );
  TRY( PetscLogEventEnd  (QPT_AllInOne,qp,0,0,0) );
  PetscFunctionReturnI(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPTFromOptions"
PetscErrorCode QPTFromOptions(QP qp)
{
  MatInvType invType=MAT_INV_MONOLITHIC;
  PetscBool ddm=PETSC_FALSE, dual=PETSC_FALSE, feti=PETSC_FALSE, project=PETSC_FALSE;
  PetscReal penalty=0.0;
  PetscBool penalty_direct=PETSC_FALSE;
  PetscBool regularize=PETSC_TRUE;

  PetscFunctionBegin;
  _fllop_ierr = PetscOptionsBegin(PetscObjectComm((PetscObject)qp),NULL,"QP transforms options","QP");CHKERRQ(_fllop_ierr);
  {
    TRY( PetscOptionsBool("-feti","perform FETI DDM combination of transforms","QPTAllInOne",feti,&feti,NULL) );
    if (feti) {
      ddm             = PETSC_TRUE;
      dual            = PETSC_TRUE;
      project         = PETSC_TRUE;
    }
    TRY( PetscOptionsReal("-penalty","QPTEnforceEqByPenalty penalty parameter","QPTEnforceEqByPenalty",penalty,&penalty,NULL) );
    TRY( PetscOptionsBool("-penalty_direct","","QPTEnforceEqByPenalty",penalty_direct,&penalty_direct,NULL) );
    TRY( PetscOptionsBool("-project","perform QPTEnforceEqByProjector","QPTEnforceEqByProjector",project,&project,NULL) );
    TRY( PetscOptionsBool("-dual","perform QPTDualize","QPTDualize",dual,&dual,NULL) );
    TRY( PetscOptionsBool("-ddm","domain decomposed data","QPTDualize",ddm,&ddm,NULL) );
    TRY( PetscOptionsBool("-regularize","perform stiffness matrix regularization (for singular TFETI matrices)","QPTDualize",regularize,&regularize,NULL) );
    if (ddm) {
      invType = MAT_INV_BLOCKDIAG;
    }
  }
  _fllop_ierr = PetscOptionsEnd();CHKERRQ(_fllop_ierr);
  TRY( QPTAllInOne(qp, invType, dual, project, penalty, penalty_direct, regularize) );
  PetscFunctionReturn(0);
}
#undef QPTransformBegin

