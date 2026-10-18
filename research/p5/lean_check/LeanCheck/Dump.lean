import Lean
open Lean Elab Command

/-! `#dump_statement thm "file"`: write the statement of `thm` and every `OAI.*` constant it depends on
(type and value, as raw kernel terms) to a file. Run once on the challenge and once on the solution, then diff. -/

namespace LeanCheck

def kindOf : ConstantInfo → String
  | .defnInfo _ => "def" | .thmInfo _ => "theorem" | .axiomInfo _ => "axiom" | .inductInfo _ => "inductive"
  | .ctorInfo _ => "ctor" | .recInfo _ => "rec" | .opaqueInfo _ => "opaque" | .quotInfo _ => "quot"

def fmt (e : Expr) : String := s!"{e.dbgToString}  #hash {e.hash}"

elab "#dump_statement " thm:ident file:str : command => do
  let env ← getEnv
  let thmName := thm.getId
  let some ci := env.find? thmName | throwError "unknown {thmName}"
  let mut entries : Array (String × String) := #[(toString thmName,
    s!"THEOREM {thmName} {ci.levelParams} : {fmt ci.type}")]
  let mut todo : Array Name := ci.type.getUsedConstants
  let mut seen : NameSet := {}
  while todo.size > 0 do
    let n := todo.back!
    todo := todo.pop
    if seen.contains n || !(`OAI).isPrefixOf n then continue
    seen := seen.insert n
    let some c := env.find? n | throwError "missing {n}"
    let mut txt := s!"{kindOf c} {n} {c.levelParams} : {fmt c.type}"
    todo := todo ++ c.type.getUsedConstants
    match c with
    | .defnInfo d => txt := txt ++ s!"\n  := {fmt d.value}"; todo := todo ++ d.value.getUsedConstants
    | .thmInfo t => txt := txt ++ "\n  := <proof omitted>"; todo := todo ++ t.value.getUsedConstants
    | .opaqueInfo o => txt := txt ++ s!"\n  := {fmt o.value}"; todo := todo ++ o.value.getUsedConstants
    | .inductInfo i => txt := txt ++ s!"\n  ctors {i.ctors}"; todo := todo ++ i.ctors.toArray
    | _ => pure ()
    entries := entries.push (toString n, txt)
  let sorted := entries.qsort (fun a b => a.1 < b.1)
  IO.FS.writeFile file.getString (String.intercalate "\n" (sorted.toList.map (·.2)) ++ "\n")
  logInfo m!"{sorted.size} constants written to {file.getString}"

end LeanCheck
